# Local pocket detail. The two overviews retain the entire ligand.
color change rgb green 0.043 0.42 0.39
color change rgb purple 0.506 0.263 0.56
color change rgb magenta 0.64 0.26 0.125
color change rgb cyan 0.22 0.48 0.65
color change rgb gray 0.075 0.29 0.48
color change rgb silver 0.58 0.65 0.72
material add PocketSurface
material change opacity PocketSurface 0.28
material change ambient PocketSurface 0.32
material change diffuse PocketSurface 0.60
material change specular PocketSurface 0.05
material add PocketRibbon
material change opacity PocketRibbon 0.90
material change ambient PocketRibbon 0.35
material change diffuse PocketRibbon 0.64
material change specular PocketRibbon 0.11

proc atlas_pocket {pdb base scale ligbonds headindices resids links {shell {}} {ribbon {}}} {
  set m [mol new $pdb type pdb waitfor all]
  set l [atomselect $m {chain B}]
  if {[$l num] != [llength $ligbonds]} {error "Ligand atom count mismatch: $pdb"}
  $l setbonds $ligbonds
  $l delete
  mol delrep 0 $m
  color change rgb silver 0.58 0.65 0.72
  if {[llength $shell] > 0} {
    mol representation QuickSurf 1.0 0.65 0.60 2.0
    mol color ColorID [colorinfo index silver]
    mol material PocketSurface
    mol selection "not chain B and index $shell"
    mol addrep $m
  }
  if {[llength $ribbon] > 0} {
    # Assign secondary structure on the intact chain, then display only the residues
    # that contribute atoms to the pocket surface.
    mol ssrecalc $m
    color change rgb cyan 0.34 0.61 0.80
    mol representation NewCartoon 0.235 18 4.0 0
    mol color ColorID [colorinfo index cyan]
    mol material PocketRibbon
    mol selection "protein and not chain B and resid $ribbon"
    mol addrep $m
  }
  if {[llength $resids] > 0} {
    if {[llength $shell] > 0} {
      mol representation Licorice 0.25 16.0 16.0
      mol color ColorID [colorinfo index gray]
    } else {
      mol representation Licorice 0.20 16.0 16.0
      mol color ColorID [colorinfo index cyan]
    }
    mol selection "not chain B and resid $resids"
    mol material ReviewStick
    mol addrep $m
  }
  mol representation Licorice 0.29 16.0 16.0
  mol color Element
  mol selection "chain B and index $headindices"
  mol material ReviewStick
  mol addrep $m

  graphics $m material ReviewStick
  foreach rgb {green purple magenta} pair $links {
    if {[llength $pair] != 2} {continue}
    graphics $m color $rgb
    lassign $pair a b
    if {[llength $shell] > 0} {
      graphics $m sphere $a radius 0.30 resolution 14
    }
    set vec [vecsub $b $a]
    set length [veclength $vec]
    if {$length < 0.0001} {continue}
    set step [expr {min(0.26, 0.70 / $length)}]
    for {set s 0.10} {$s < 0.92} {set s [expr {$s+$step}]} {
      set f [expr {min(0.92,$s+0.59*$step)}]
      set radius [expr {[llength $shell] > 0 ? 0.14 : 0.105}]
      graphics $m cylinder [vecadd $a [vecscale $s $vec]] [vecadd $a [vecscale $f $vec]] radius $radius resolution 10 filled yes
    }
  }
  display resetview
  set zoom [list [list $scale 0 0 0] [list 0 $scale 0 0] [list 0 0 $scale 0] {0 0 0 1}]
  molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} [list [transidentity] [transidentity] $zoom [transidentity]]
  render Tachyon $base.dat true
  if {![file exists $base.dat] || [file size $base.dat] < 1000} {error "Missing Tachyon scene $base"}
  mol delete $m
  puts "ATLAS_RENDERED [file tail $base]"
}
