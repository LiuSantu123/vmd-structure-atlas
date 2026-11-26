"""Compose paginated A3 cards using the VMD views and measured contacts."""
from collections import defaultdict
import io
import math

import numpy as np
from PIL import Image
import pymupdf as pdf

W, H = 297*72/25.4, 420*72/25.4
INK, SOFT, BG, LINE = (.06, .15, .23), (.27, .35, .43), (.95, .97, .98), (.85, .89, .93)
LINKS = ((.043, .42, .39), (.506, .263, .56), (.64, .26, .125))
NORMAL, BOLD = 'helv', 'hebo'
TYPES = {'basic–phosphate': 'BP', 'polar proximity': 'POL', 'C–C proximity': 'CC',
         'other proximity': 'NEAR', 'short overlap': 'SHORT'}


def text(page, x, y, s, size=8, color=INK, bold=False, bounds=None, right=False):
    font = BOLD if bold else NORMAL
    s = str(s).encode('cp1252', errors='replace').decode('cp1252')
    limit = bounds if bounds is not None else W
    width = pdf.get_text_length(s, fontname=font, fontsize=size)
    if width > limit:
        # Shorten only the displayed label; full text is in the index CSV.
        shortened = s
        while shortened and pdf.get_text_length(shortened+'…', fontname=font, fontsize=size) > limit:
            shortened = shortened[:-1]
        s = shortened+'…'
        width = pdf.get_text_length(s, fontname=font, fontsize=size)
    page.insert_text((x-width if right else x, y), s, fontname=font, fontsize=size, color=color)


def place(page, image, crop, frame):
    l, t, r, b = crop
    x, y, w, h = frame
    f = min(w/(r-l), h/(b-t))
    ox, oy = x+(w-f*(r-l))/2, y+(h-f*(b-t))/2
    stream = io.BytesIO()
    image.crop(crop).save(stream, 'JPEG', quality=92, subsampling=0)
    page.insert_image(pdf.Rect(ox, oy, ox+f*(r-l), oy+f*(b-t)), stream=stream.getvalue())
    return lambda q: (ox+(q[0]-l)*f, oy+(q[1]-t)*f)


def framing(cards):
    """One crop for both overview angles and all members of an aligned group."""
    boxes = defaultdict(list)
    for card in cards:
        boxes[card['group']].extend(card['boxes'][:2])
    return {group: (max(0, min(b[0] for b in bs)-12), max(0, min(b[1] for b in bs)-12),
                    min(654, max(b[2] for b in bs)+12), min(422, max(b[3] for b in bs)+12))
            for group, bs in boxes.items()}


def badges(image, shown):
    occupied = np.any(np.asarray(image.convert('RGB')) < 220, axis=-1)
    h, w = occupied.shape
    points = []
    for a in shown:
        # Reserve edge space even when the dashed line lies near the image edge.
        x, y = np.clip(a['mid_pixel'], (18, 18), (w-19, h-19))
        choices = []
        for d in (28, 38, 50, 64):
            for t in np.arange(0, 2*np.pi, np.pi/8):
                xx, yy = round(x+d*np.cos(t)), round(y+d*np.sin(t))
                r = 13
                if r < xx < w-r and r < yy < h-r:
                    ink = occupied[yy-r:yy+r+1, xx-r:xx+r+1].mean()
                    overlap = sum(max(0, 32-math.hypot(xx-x2, yy-y2)) for x2, y2 in points)
                    choices.append((95*ink + .2*(d-28) + 2*overlap, xx, yy))
        points.append(min(choices)[1:] if choices else (int(x), int(y)))
    return points


def rule(page, x, y, xx, yy):
    page.draw_line((x, y), (xx, yy), color=LINE, width=.45)


def card_view(page, card, root, x, y, cw, pitch, crop, cutoff):
    header = 29
    ref = str(card['id'])
    if card['label']:
        ref += '  /  '+card['label']
    elif card['ligand_code']:
        ref += '  /  '+card['ligand_code']
    meta = '  ·  '.join(s for s in (card['model'], f"seed {card['seed']}" if card['seed'] else '') if s) or 'model / seed unrecorded'
    head = card['head_residue_count']
    other = card['other_residue_count']
    meta += f"   H {head}" + (f' / O {other}' if other is not None else '') + (' !' if card['n_short_pairs'] else '')
    page.draw_rect(pdf.Rect(x, y, x+cw, y+header), fill=BG, color=None)
    text(page, x+7, y+12, ref, 10.15, bold=True, bounds=cw*.58)
    text(page, x+cw-7, y+12, meta, 8.15, SOFT, bounds=cw*.37, right=True)
    metrics = [('iptm', 'iPTM', 2, ''), ('ipsae', 'ipSAE', 3, ''), ('mean_pae', 'Mean PAE', 1, ' Å'),
               ('protein_plddt', 'Prot pLDDT', 1, ''), ('ligand_plddt', 'Lig pLDDT', 1, '')]
    for frac, (key, label, n, unit) in zip((.018, .18, .36, .58, .81), metrics):
        v = card['values'][key]
        s = '—' if v is None else f'{v:.{n}f}{unit}'
        text(page, x+cw*frac, y+23, f'{label} {s}', 7.95, SOFT)
    end_y = y+pitch-16
    top = y+header+1
    edges = (0, .285*cw, .57*cw, cw)
    for j, name in enumerate(card['images']):
        if j:
            rule(page, x+edges[j], top+3, x+edges[j], end_y-3)
        frame = (x+edges[j]+1, top, edges[j+1]-edges[j]-2, end_y-top-2)
        with Image.open(root/'views'/name) as im:
            if j < 2:
                place(page, im, crop, frame)
                continue
            positions = badges(im, card['shown'])
            l, t, r, b = card['boxes'][2]
            pocket_crop = (max(0, min([l-12]+[xx-16 for xx, _ in positions])),
                           max(0, min([t-12]+[yy-16 for _, yy in positions])),
                           min(im.width, max([r+12]+[xx+16 for xx, _ in positions])),
                           min(im.height, max([b+12]+[yy+16 for _, yy in positions])))
            convert = place(page, im, pocket_crop, frame)
        for a, label in zip(card['shown'], positions):
            ix = a['badge']
            c = LINKS[ix-1]
            xx, yy = convert(label)
            mid = convert(a['mid_pixel'])
            page.draw_line(mid, (xx, yy), color=c, width=.50)
            page.draw_circle((xx, yy), radius=4.1, color=c, fill=(1, 1, 1), width=.60)
            text(page, xx-1.65, yy+2.15, ix, 6.2, c, bold=True)
    rule(page, x+3, end_y, x+cw-3, end_y)
    for i, a in enumerate(card['shown']):
        cid = f"{a['protein_chain']}:" if a['protein_chain'] else ''
        s = (f"{a['badge']}  {cid}{a['protein_resname']}{a['protein_resid']} {a['protein_atom']}–"
             f"{a['ligand_atom']}  {a['distance_A']:.1f} Å  {TYPES[a['class']]}")
        text(page, x+6+i*(cw-9)/3, end_y+10, s, 7.1, LINKS[i], bold=True, bounds=(cw-12)/3)
    if not card['shown']:
        text(page, x+6, end_y+10, f'No protein–head heavy atom pair within {cutoff:g} Å', 7.3, SOFT)


def compose_pdf(cards, output, opts):
    doc = pdf.open()
    per_page = opts.rows_per_page*2
    left, gap, card_width = 13, 10, (W-36)/2
    pitch = min((H-152)/opts.rows_per_page, 365)
    crops = framing(cards)
    bookmark = []
    for start in range(0, len(cards), per_page):
        items = cards[start:start+per_page]
        page = doc.new_page(width=W, height=H)
        number = page.number+1
        text(page, left, 34, opts.title, 21, bold=True, bounds=W-130)
        text(page, W-left, 34, f'PAGE {number}', 9, SOFT, bold=True, right=True)
        text(page, left, 53, f'Cards {start+1}–{start+len(items)} of {len(cards)}  ·  VMD / Tachyon', 9, SOFT)
        rule(page, left, 60, W-left, 60)
        bookmark.append([1, f"{items[0]['id']} — {items[-1]['id']}" if len(items)>1 else items[0]['id'], number])
        for c in range(min(len(items), 2)):
            x = left+c*(card_width+gap)
            text(page, x+3, 72, 'A  BINDING FACE', 8, SOFT, True)
            text(page, x+.285*card_width+3, 72, 'B  SIDE / 82°', 8, SOFT, True)
            text(page, x+.57*card_width+3, 72, 'C  POCKET SURFACE + CARTOON', 8, SOFT, True)
        for i, card in enumerate(items):
            r, c = divmod(i, 2)
            card.update(pdf_page=number, pdf_card=i+1)
            card_view(page, card, output, left+c*(card_width+gap), 79+r*pitch, card_width, pitch, crops[card['group']], opts.contact)
        rule(page, left, H-63, W-left, H-63)
        text(page, left, H-49, f'A/B: full pose and 75% transparent surface. C: {opts.pocket_radius:g} Å pocket shell and corresponding ribbon; sticks mark selected residues.', 8, SOFT)
        text(page, left, H-37, f'H/O: residue counts within {opts.contact:g} Å of head/other ligand atoms. BP: basic–phosphate; POL: N/O; CC: C–C; NEAR: other; !/SHORT: pair <2.2 Å.', 8, SOFT)
        text(page, left, H-25, 'Links show geometric proximity. —: unavailable. Metric sources and complete contact list are in the CSV files.', 8, SOFT)
    doc.set_toc(bookmark)
    doc.set_metadata({'title': opts.title, 'subject': 'Protein–ligand overviews, pocket cartoon and geometric contacts'})
    doc.save(output/'.atlas.partial.pdf', garbage=3, deflate=True)
    doc.close()
    (output/'.atlas.partial.pdf').replace(output/'atlas.pdf')
