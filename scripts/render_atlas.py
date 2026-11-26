#!/usr/bin/env python3
"""Render a CSV of existing protein–ligand structures with VMD and Tachyon."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import gemmi
import numpy as np
from PIL import Image, ImageChops
from scipy.spatial import cKDTree, distance

from compose import compose_pdf

HERE = Path(__file__).resolve().parent
RAY = (654, 422)
METRICS = ('iptm', 'ipsae', 'mean_pae', 'protein_plddt', 'ligand_plddt')
CONTACT_FIELDS = ('id', 'protein_chain', 'protein_resid', 'protein_resname', 'protein_atom',
                  'ligand_atom', 'region', 'distance_A', 'class', 'badge')
PROTEIN_CHAINS = 'ACDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
BASIC = {'ARG': {'NE', 'NH1', 'NH2'}, 'LYS': {'NZ'}}


@dataclass
class Pose:
    atoms: list
    ligand: list
    residues: list
    atom_resid: np.ndarray  # zero-based index into residues
    alpha: np.ndarray      # indices of the Cα atoms in atoms
    ligand_chain: str
    ligand_code: str
    ligand_resid: str
    xyz: np.ndarray

    @property
    def ca(self):
        return self.xyz[self.alpha]

    @property
    def signature(self):
        return tuple((c, str(r.seqid), r.name) for c, r in self.residues)


def heavy_atoms(residue):
    # Blank altloc is shared; otherwise choose one conformer for the whole residue.
    others = defaultdict(float)
    for a in residue:
        if a.altloc not in ('\x00', ' '):
            others[a.altloc] += a.occ
    alt = min(others, key=lambda k: (-others[k], k != 'A', k)) if others else None
    atoms = {}
    for a in residue:
        if a.element.name not in ('H', 'D') and a.altloc in ('\x00', ' ', alt):
            if a.name not in atoms or a.altloc in ('\x00', ' '):
                atoms[a.name] = a
    return list(atoms.values())


def load_pose(row):
    path = Path(row['structure'])
    model = gemmi.read_structure(str(path))[0]
    chains = set(row['protein_chains'].split())
    atoms, residues, ridx, alpha, options = [], [], [], [], []
    for chain in model:
        for res in chain:
            aa = gemmi.find_tabulated_residue(res.name).is_amino_acid()
            heavy = heavy_atoms(res)
            if aa:
                if chains and chain.name not in chains:
                    continue
                ca = [i for i, a in enumerate(heavy) if a.name == 'CA']
                if not ca:
                    continue
                alpha.append(len(atoms) + ca[0])
                ridx.extend([len(residues)] * len(heavy))
                atoms.extend(heavy)
                residues.append((chain.name, res))
            elif heavy and not res.is_water() and all(not row[k] or row[k] == v for k, v in (
                    ('ligand_code', res.name), ('ligand_chain', chain.name), ('ligand_resid', str(res.seqid)))):
                options.append((chain.name, res, heavy))
    if len(options) != 1:
        keys = ', '.join(f'{c or "(blank)"}:{r.seqid}:{r.name}' for c, r, _ in options[:12])
        raise ValueError(f"Expected one ligand residue, found {len(options)} ({keys}); specify ligand_code/chain/resid")
    if chains - {c for c, _ in residues} or len(residues) < 3:
        raise ValueError('Protein selection must include at least three Cα residues and all named chains')
    c, r, lig = options[0]
    if len(residues) > 9999 or len(atoms) + len(lig) > 99999 or len(lig) > 999 or len({c for c, _ in residues}) > len(PROTEIN_CHAINS):
        raise ValueError('Complex exceeds PDB staging limits; select fewer protein chains')
    xyz = np.array([a.pos.tolist() for a in atoms + lig], dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError('Non-finite coordinates')
    return Pose(atoms, lig, residues, np.array(ridx), np.array(alpha), c, r.name, str(r.seqid), xyz)


def bond_pairs(path, code, explicit=False):
    doc = gemmi.cif.read_file(str(path))
    output = []
    for block in doc:
        table = block.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2', '?comp_id'])
        for entry in table:
            comp = gemmi.cif.as_string(entry[2]) if entry.has(2) else ''
            if comp == code or (explicit and not comp and len(doc) == 1):
                output.append(tuple(gemmi.cif.as_string(entry[i]) for i in (0, 1)))
    return output


def connected(bonds, n):
    neighbors = [set() for _ in range(n)]
    for i, j in bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
    seen, queue = set(), [0]
    while queue:
        i = queue.pop()
        if i not in seen:
            seen.add(i)
            queue.extend(neighbors[i] - seen)
    return len(seen) == n


def ligand_geometry(row, pose):
    names = {a.name: i for i, a in enumerate(pose.ligand)}
    head_names = row['head_atoms'].split()
    if head_names and not set(head_names) <= names.keys():
        raise ValueError(f'Unknown ligand head atom(s): {sorted(set(head_names) - names.keys())}')
    head = np.array(sorted({names[a] for a in head_names} if head_names else names.values()), dtype=int)
    if len(names) == 1:
        return head, [], 'single atom'
    source = Path(row['structure'])
    candidates = []
    if row['ccd']:
        candidates.append(('supplied CCD', bond_pairs(row['ccd'], pose.ligand_code, explicit=True)))
    elif source.suffix.lower() in ('.cif', '.mmcif'):
        candidates.append(('embedded CCD', bond_pairs(source, pose.ligand_code)))
    else:
        serial = {a.serial: a.name for a in pose.ligand}
        pairs = []
        with source.open(errors='replace') as f:
            for line in f:
                if line.startswith('CONECT'):
                    fields = [line[i:i+5].strip() for i in range(6, len(line.rstrip()), 5)]
                    ids = [int(x) for x in fields if x.isdigit()]
                    if ids and ids[0] in serial:
                        pairs.extend((serial[ids[0]], serial[j]) for j in ids[1:] if j in serial)
        candidates.append(('PDB CONECT', pairs))
    for provenance, pairs in candidates:
        bonds = sorted({tuple(sorted((names[a], names[b]))) for a, b in pairs if a != b and a in names and b in names})
        if connected(bonds, len(names)):
            return head, bonds, provenance
    if row['ccd']:
        raise ValueError('Supplied CCD does not define a connected graph for these ligand atom names')
    if row['infer_bonds'].lower() in ('1', 'true', 'yes'):
        rad = np.array([a.element.covalent_r for a in pose.ligand])
        d = distance.squareform(distance.pdist(pose.xyz[len(pose.atoms):]))
        keep = np.triu((d > 0.65) & (d < 1.23 * (rad[:, None] + rad)), k=1)
        bonds = [tuple(map(int, ij)) for ij in np.argwhere(keep)]
        if connected(bonds, len(names)):
            return head, bonds, 'inferred from length; inspect'
    raise ValueError('Ligand graph unavailable or disconnected; supply a matching CCD or explicitly set infer_bonds=true')


def confidence(row, pose):
    values = {k: float(row[k]) if row[k] else None for k in METRICS}
    notes = {k: 'manifest' for k, v in values.items() if v is not None}
    if row['af3_confidence']:
        data = json.loads(Path(row['af3_confidence']).read_bytes())
        pc = {c for c, _ in pose.residues}
        lc = pose.ligand_chain
        if lc in pc:
            raise ValueError('AF3 confidence requires distinct protein and ligand chains; provide precomputed values for a shared chain')
        ids = np.asarray(data['token_chain_ids'])
        p, l = np.isin(ids, list(pc)), ids == lc
        pae = np.asarray(data['pae'], dtype=float)
        if (pae.shape != (len(ids), len(ids)) or not np.isfinite(pae).all()
                or np.any(pae < 0) or any(np.sum(ids == c) != sum(a == c for a, _ in pose.residues) for c in pc)
                or p.sum() != len(pose.residues) or l.sum() != len(pose.ligand)):
            raise ValueError('AF3 PAE has no matching protein-residue / ligand-heavy-atom token mapping')
        blocks = (pae[np.ix_(p, l)], pae[np.ix_(l, p)])

        def ipsae(block):
            # Dunbrack d0res, applied to AF3 ligand heavy-atom tokens.
            mask = block < 10.0
            n = mask.sum(axis=1)
            d0 = np.maximum(1.0, 1.24 * np.cbrt(np.maximum(n, 26) - 15) - 1.8)
            numer = (mask / (1.0 + (block / d0[:, None])**2)).sum(axis=1)
            return float(np.max(np.divide(numer, n, out=np.zeros(len(n)), where=n > 0)))

        found = {'mean_pae': float((blocks[0].mean() + blocks[1].mean()) / 2),
                 'ipsae': max(ipsae(x) for x in blocks)}
        if 'atom_plddts' in data:
            ac = np.asarray(data['atom_chain_ids'])
            pl = np.asarray(data['atom_plddts'], dtype=float)
            if pl.shape != ac.shape or not np.isfinite(pl).all() or np.any((pl < 0) | (pl > 100.001)):
                raise ValueError('Invalid AF3 atom pLDDT')
            for k, mask in [('protein_plddt', np.isin(ac, list(pc))), ('ligand_plddt', ac == lc)]:
                if mask.any():
                    found[k] = float(pl[mask].mean())
        for k, v in found.items():
            if values[k] is None:
                values[k] = v
                notes[k] = 'AF3 d0res max direction, PAE < 10 A' if k == 'ipsae' else 'AF3 full confidence'
    if row['af3_summary'] and values['iptm'] is None:
        v = json.loads(Path(row['af3_summary']).read_bytes()).get('iptm')
        if v is not None:
            values['iptm'], notes['iptm'] = float(v), 'AF3 summary (overall iPTM)'
    for k, v in values.items():
        ceiling = 100.001 if k.endswith('plddt') else 1.001 if k in ('iptm', 'ipsae') else np.inf
        if v is not None and not (np.isfinite(v) and 0 <= v <= ceiling):
            raise ValueError(f'Invalid {k}: {v}; use a blank cell for missing data')
    return values, '; '.join(f'{k}: {v}' for k, v in notes.items())


def unit(x):
    return x / np.linalg.norm(x)


def basis(ca, outward):
    _, axes = np.linalg.eigh((ca - ca.mean(axis=0)).T @ (ca - ca.mean(axis=0)))
    z = unit(outward) if np.linalg.norm(outward) > 1e-6 else axes[:, 0]
    y = next(x for x in axes.T[::-1] if np.linalg.norm(np.cross(x, z)) > 0.25)
    y = unit(y - z * (y @ z))
    if y @ (ca[-1] - ca[0]) < 0:
        y = -y
    return np.column_stack((np.cross(y, z), y, z))


def ry(deg):
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.array(((c, 0, s), (0, 1, 0), (-s, 0, c)))


def rx(deg):
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.array(((1, 0, 0), (0, c, -s), (0, s, c)))


def kabsch(mobile, fixed):
    u, _, vt = np.linalg.svd((mobile - mobile.mean(axis=0)).T @ (fixed - fixed.mean(axis=0)))
    return u @ np.diag((1., 1., np.linalg.det(u @ vt))) @ vt


def overview_plan(rows):
    states, result = {}, {}
    spans = defaultdict(float)
    for row in rows:
        try:
            pose = load_pose(row)
            names = set(row['head_atoms'].split())
            ligand_xyz = pose.xyz[len(pose.atoms):]
            head = ligand_xyz[[a.name in names for a in pose.ligand]] if names else ligand_xyz
            if names - {a.name for a in pose.ligand}:
                raise ValueError('head_atoms contains names absent from the ligand')
            group = ('group', row['align_group']) if row['align_group'] else ('card', row['id'])
            if group not in states:
                states[group] = pose.ca, pose.signature, basis(pose.ca, head.mean(axis=0) - pose.ca.mean(axis=0))
            ca, signature, rot = states[group]
            if pose.signature != signature:
                raise ValueError(f"align_group {row['align_group']} needs the same ordered protein chain/residue IDs and names")
            align = kabsch(pose.ca, ca) if not np.array_equal(pose.ca, ca) else np.eye(3)
            rots = (align @ rot, align @ rot @ ry(82))
            for tr in rots:
                half = np.ptp(pose.xyz @ tr, axis=0) / 2
                spans[group] = max(spans[group], half[0]/1.38, half[1])
            result[row['id']] = (group, rots)
        except Exception as err:
            raise ValueError(f"{row['id']}: {err}") from err
    return {i: (g, rots, 1.39 / max(16., spans[g]+3.4)) for i, (g, rots) in result.items()}


def contact_geometry(row, pose, head, bonds, cutoff, radius):
    n = len(pose.atoms)
    protein, ligand = pose.xyz[:n], pose.xyz[n:]
    is_head = set(head.tolist())
    neighbors = [set() for _ in pose.ligand]
    for i, j in bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
    pairs = []
    rank = ('basic–phosphate', 'polar proximity', 'C–C proximity', 'other proximity', 'short overlap')
    for i, js in enumerate(cKDTree(protein).query_ball_tree(cKDTree(ligand), cutoff)):
        for j in js:
            r = int(pose.atom_resid[i])
            chain, res = pose.residues[r]
            a, b = pose.atoms[i], pose.ligand[j]
            d = float(np.linalg.norm(protein[i] - ligand[j]))
            terminal_phosphate = b.element.name == 'O' and len(neighbors[j]) == 1 and any(
                pose.ligand[k].element.name == 'P' for k in neighbors[j])
            kind = (4 if d < 2.2 else 0 if d <= 3.6 and terminal_phosphate and a.name in BASIC.get(res.name, set())
                    else 1 if d <= 3.5 and {a.element.name, b.element.name} <= {'N', 'O'}
                    else 2 if a.element.name == b.element.name == 'C' else 3)
            pairs.append({'id': row['id'], 'protein_chain': chain, 'protein_resid': str(res.seqid),
                          'protein_resname': res.name, 'protein_atom': a.name, 'ligand_atom': b.name,
                          'region': 'head' if j in is_head else 'other', 'distance_A': round(d, 3),
                          'class': rank[kind], 'badge': '', 'r': r, 'p': i, 'l': j, 'rank': kind})
    pairs.sort(key=lambda a: (a['rank'], a['distance_A'], a['r'], a['p'], a['l']))
    shown, used = [], set()
    for a in pairs:
        if a['region'] == 'head' and a['r'] not in used:
            shown.append(a)
            used.add(a['r'])
            a['badge'] = len(shown)
            if len(shown) == 3:
                break
    shell = np.flatnonzero(cKDTree(ligand[head]).query(protein)[0] <= radius)
    return pairs, shown, shell


def pocket_camera(pose, head, shown, shell):
    xyz, ca = pose.xyz, pose.ca
    head_idx = head + len(pose.atoms)
    center = xyz[head_idx].mean(axis=0)
    outward = center - ca.mean(axis=0)
    if len(shell):
        local = center - xyz[shell].mean(axis=0)
        if np.linalg.norm(local) > .5 and local @ outward > 0:
            outward = .55 * unit(local) + .45 * unit(outward)
    base = basis(ca, outward)
    selected = np.array([a['p'] for a in shown], dtype=int)
    ligand = len(pose.atoms) + np.array([a['l'] for a in shown], dtype=int)
    # Frame the full annotated sticks plus Cα of the residues shown as ribbon.
    sticks = np.flatnonzero(np.isin(pose.atom_resid, [a['r'] for a in shown]))
    shell_resids = np.unique(pose.atom_resid[shell])
    fits = np.unique(np.r_[head_idx, sticks, shell, pose.alpha[shell_resids]])
    best = None
    for tilt in (-35, -17, 0, 17, 35):
        for turn in (-50, -30, -15, 0, 15, 30, 50):
            rot = base @ ry(turn) @ rx(tilt)
            z = xyz @ rot
            pa, la = z[selected, :2], z[ligand, :2]
            link = np.linalg.norm(pa - la, axis=1)
            separation = distance.pdist((pa + la) / 2)
            extent = np.ptp(z[fits, :2], axis=0)
            score = (np.minimum(link, 3.3).sum() + .58 * np.minimum(separation, 5.).sum()
                     - .068 * (extent[1] + extent[0]/1.55) - .011 * (abs(tilt)+abs(turn)))
            if len(shell):
                occlusion = (distance.cdist(z[shell, :2], z[head_idx, :2]) < 1.65) & (z[shell, 2, None] > z[head_idx, 2] + .7)
                score -= 8 * occlusion.any(axis=0).mean()
            if best is None or score > best[0]:
                best = score, z
    z = best[1]
    dx, dy = np.ptp(z[fits, :2], axis=0)
    if max(dy/1.55, dx) + .3 < max(dx/1.55, dy):
        z = z @ np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
    z -= (z[fits].max(axis=0) + z[fits].min(axis=0)) / 2
    span = np.ptp(z[fits, :2], axis=0) / 2
    scale = .83 * 1.32 / max(9.0, (span[0]+2.4)/1.55, span[1]+2.4)
    return z, scale, shell_resids


def write_pdb(path, pose, xyz):
    if xyz.min() < -999.999 or xyz.max() > 9999.999:
        raise ValueError('Coordinates exceed fixed-width PDB fields after centering')
    chains = dict(zip(dict.fromkeys(c for c, _ in pose.residues), PROTEIN_CHAINS))

    def line(i, a, name, rname, chain, rnum, record):
        x, y, z = xyz[i]
        return (f'{record:<6}{i+1:5d} {name:4s} {rname:>3s} {chain}{rnum:4d}    '
                f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.element.name:>2s}\n')

    with path.open('w') as f:
        prev = None
        for i, (a, r) in enumerate(zip(pose.atoms, pose.atom_resid)):
            chain, res = pose.residues[r]
            if prev is not None and prev != chain:
                f.write('TER\n')
            prev = chain
            name = a.name.ljust(4) if len(a.name) == 4 else ' ' + a.name.ljust(3)
            f.write(line(i, a, name, res.name, chains[chain], r+1, 'ATOM'))
        f.write('TER\n')
        for j, a in enumerate(pose.ligand):
            f.write(line(len(pose.atoms)+j, a, f'{a.element.name[:1]}{j+1:03d}', 'LIG', 'B', 1, 'HETATM'))
        f.write('END\n')


def tcl_list(values):
    return '{' + ' '.join(map(str, values)) + '}'


def trace(tachyon, work, image_dir, tag, dest):
    ppm = work / f'{tag}.ppm'
    log = work / f'{tag}.log'
    with log.open('w') as f:
        result = subprocess.run([tachyon, f'{tag}.dat', '-numthreads', '1', '-res', *map(str, RAY),
                                 '-aasamples', '12', '-trans_vmd', '-format', 'PPM', '-o', ppm.name],
                                cwd=work, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode or not ppm.exists():
        raise RuntimeError(log.read_text(errors='replace')[-3500:])
    with Image.open(ppm) as raw:
        im = raw.convert('RGB')
        box = ImageChops.difference(im, Image.new('RGB', im.size, 'white')).point(lambda p: 255 if p > 9 else 0).getbbox()
        if box is None or min(box[0], box[1], im.width-box[2], im.height-box[3]) < 3:
            raise RuntimeError(f'{tag}: blank or clipped rendered frame ({box})')
        im.save(image_dir / dest, 'PNG')
    return box


def render_card(row, plan, opts, vmd, tachyon, output, number):
    try:
        pose = load_pose(row)
        head, bonds, topology = ligand_geometry(row, pose)
        values, notes = confidence(row, pose)
        n = len(pose.atoms)
        pairs, shown, shell = contact_geometry(row, pose, head, bonds, opts.contact, opts.pocket_radius)
        pocket, zoom, cartoon = pocket_camera(pose, head, shown, shell)
        adjacency = [[] for _ in pose.ligand]
        for a, b in bonds:
            adjacency[a].append(n+b)
            adjacency[b].append(n+a)
        neighbors = tcl_list(tcl_list(a) for a in adjacency)
        group, rotations, wide_scale = plan
        tags = ('binding-face', 'side', 'pocket')
        names = tuple(f'{number:05d}_{t}.png' for t in tags)
        with tempfile.TemporaryDirectory(prefix='vmd-atlas-') as tmp:
            work = Path(tmp)
            commands = []
            for tag, rot in zip(tags[:2], rotations):
                z = pose.xyz @ rot
                z -= (z.max(axis=0) + z.min(axis=0))/2
                write_pdb(work/f'{tag}.pdb', pose, z)
                commands.append(f'atlas_overview {tag}.pdb {tag} {wide_scale:.9f} {neighbors}')
            write_pdb(work/'pocket.pdb', pose, pocket)
            links = []
            for a in shown:
                z1, z2 = pocket[a['p']], pocket[n+a['l']]
                links.append(tcl_list(tcl_list(x) for x in (z1, z2)))
                mid = (z1+z2)/2
                a['mid_pixel'] = (np.array(RAY)/2 + mid[:2] * np.array((1, -1)) * zoom*RAY[1]/3).tolist()
            commands.append(' '.join(('atlas_pocket pocket.pdb pocket', f'{zoom:.9f}', neighbors,
                                      tcl_list(n + head), tcl_list(a['r']+1 for a in shown), tcl_list(links),
                                      tcl_list(shell), tcl_list(cartoon+1))))
            (work/'jobs.tcl').write_text('\n'.join(commands) + '\n')
            (work/'launch.tcl').write_text('if {[catch {\nsource $::env(VMD_ATLAS_OVERVIEW)\n'
                                           'source $::env(VMD_ATLAS_POCKET)\nsource jobs.tcl\n'
                                           '} reason]} {puts stderr $errorInfo; exit 1}\nquit\n')
            env = dict(os.environ, VMDNOCUDA='1', VMDNOOPTIX='1', VMDFORCECPUCOUNT='2',
                       VMD_ATLAS_OVERVIEW=str(HERE/'overview.tcl'), VMD_ATLAS_POCKET=str(HERE/'pocket.tcl'))
            log = output/'logs'/f'{number:05d}.log'
            with log.open('w') as f:
                done = subprocess.run([vmd, '-dispdev', 'text', '-e', 'launch.tcl'], cwd=work, env=env,
                                      stdin=subprocess.DEVNULL, stdout=f, stderr=subprocess.STDOUT)
            text = log.read_text(errors='replace')
            if done.returncode != 0 or len(re.findall(r'^ATLAS_RENDERED ', text, flags=re.M)) != len(tags):
                raise RuntimeError(f'VMD did not produce all three scenes; see {log}\n{text[-2600:]}')
            boxes = [trace(tachyon, work, output/'views', tag, name) for tag, name in zip(tags, names)]
        return dict(row, **values) | dict(number=number, values=values, metric_sources=notes, topology=topology,
                    ligand_selection=f'{pose.ligand_chain}:{pose.ligand_resid}:{pose.ligand_code}',
                    images=names, boxes=boxes, group=group, shown=shown, contacts=pairs,
                    shell_residues=';'.join(f'{pose.residues[i][0]}:{pose.residues[i][1].seqid}' for i in cartoon),
                    head_residue_count=len({a['r'] for a in pairs if a['region'] == 'head'}),
                    other_residue_count=len({a['r'] for a in pairs if a['region'] == 'other'}) if len(head)<len(pose.ligand) else None,
                    n_short_pairs=sum(a['rank'] == 4 for a in pairs),
                    structure_sha256=hashlib.sha256(Path(row['structure']).read_bytes()).hexdigest())
    except Exception as err:
        raise RuntimeError(f"{row['id']} ({row['structure']}): {err}") from err


def resolve_executable(value, candidates, label):
    for name in [value] if value else candidates:
        if name and (p := shutil.which(str(name))):
            return os.path.abspath(p)
    raise SystemExit(f'{label} executable not found; pass --{label.lower()} or use an installed VMD prefix')


def read_manifest(path):
    required = ('id', 'structure')
    optional = ('label', 'model', 'seed', 'ligand_code', 'ligand_chain', 'ligand_resid', 'protein_chains',
                'head_atoms', 'ccd', 'infer_bonds', 'align_group', 'af3_confidence', 'af3_summary', *METRICS)
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError('Manifest requires the id and structure columns')
        rows = []
        for line, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f'Manifest line {line}: unexpected extra cells')
            row = {k: (raw.get(k) or '').strip() for k in (*required, *optional)}
            if not all(row[k] for k in required):
                raise ValueError(f'Manifest line {line}: id and structure must not be blank')
            for k in ('structure', 'ccd', 'af3_confidence', 'af3_summary'):
                if row[k]:
                    p = (path.parent/Path(row[k]).expanduser()).resolve()
                    if not p.is_file():
                        raise ValueError(f'Manifest line {line}: no {k} file at {p}')
                    row[k] = str(p)
            rows.append(row)
    if not rows or len(set(r['id'] for r in rows)) != len(rows):
        raise ValueError('Manifest must have at least one row and unique card IDs')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('--out', required=True, type=Path, help='new output directory')
    p.add_argument('--title', default='Protein–ligand structure atlas')
    p.add_argument('--preview', type=int, metavar='N', help='render only the first N manifest rows')
    p.add_argument('--rows-per-page', type=int, default=8, choices=range(1, 9), metavar='1–8')
    p.add_argument('--contact', type=float, default=4., help='heavy atom cutoff in Å')
    p.add_argument('--pocket-radius', type=float, default=8.5, help='pocket shell radius in Å')
    p.add_argument('--workers', type=int, default=2, help='independent VMD/Tachyon card workers')
    p.add_argument('--vmd')
    p.add_argument('--tachyon')
    opts = p.parse_args()
    if opts.workers < 1 or (opts.preview is not None and opts.preview < 1) or not (0 < opts.contact <= opts.pocket_radius < np.inf):
        p.error('Require positive workers/preview and 0 < contact <= pocket-radius')
    vmd = resolve_executable(opts.vmd or os.environ.get('VMD_ATLAS_VMD'), ['vmd'], 'VMD')
    prefix = Path(vmd).parent.parent
    bundled_dirs = (prefix/'lib'/'vmd', Path(vmd).parent, Path(os.environ.get('VMDDIR', '/nonexistent')))
    tachyon = resolve_executable(opts.tachyon or os.environ.get('VMD_ATLAS_TACHYON'),
                                 ['tachyon', *(p for d in bundled_dirs for p in sorted(d.glob('tachyon_*')))], 'Tachyon')
    rows = read_manifest(opts.manifest.resolve())[:opts.preview]
    plan = overview_plan(rows)
    output = opts.out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output/'logs').mkdir()
    (output/'views').mkdir()

    def render(item):
        number, row = item
        return render_card(row, plan[row['id']], opts, vmd, tachyon, output, number)

    cards = []
    with ThreadPoolExecutor(max_workers=opts.workers) as workers:
        for card in workers.map(render, enumerate(rows, 1)):
            cards.append(card)
            print(f"{len(cards)}/{len(rows)}  {card['id']}", flush=True)
    compose_pdf(cards, output, opts)
    fields = ('id', 'label', 'model', 'seed', 'align_group', 'pdf_page', 'pdf_card', 'structure',
              'structure_sha256', 'ligand_selection', 'head_atoms', 'topology', *METRICS, 'metric_sources',
              'head_residue_count', 'other_residue_count', 'n_short_pairs', 'shell_residues')
    for filename, keys, records in (('structure_index.csv', fields, cards),
                                    ('heavy_atom_contacts.csv', CONTACT_FIELDS, (a for c in cards for a in c['contacts']))):
        with (output/filename).open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(records)
    print(output/'atlas.pdf', flush=True)


if __name__ == '__main__':
    main()
