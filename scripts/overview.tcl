# Build VMD scenes for actual CPU ray tracing with Tachyon.
display projection Orthographic
display depthcue off
display culling off
display shadows off
display ambientocclusion off
display antialias on
display resize 530 342
axes location Off
color Display Background white
color change rgb blue 0.14 0.39 0.64
color change rgb silver 0.63 0.76 0.87
color change rgb orange 0.73 0.31 0.07
color change rgb yellow 0.98 0.66 0.04
color change rgb red 0.83 0.08 0.11
color Element C orange
color Element O red
color Element P yellow
color Element N blue
material add ReviewEnvelope
material change opacity ReviewEnvelope 0.25
material change ambient ReviewEnvelope 0.32
material change diffuse ReviewEnvelope 0.63
material change specular ReviewEnvelope 0.10
material add ReviewCartoon
material change opacity ReviewCartoon 1.00
material change ambient ReviewCartoon 0.24
material change diffuse ReviewCartoon 0.74
material change specular ReviewCartoon 0.18
material add ReviewStick
material change opacity ReviewStick 1.00
material change ambient ReviewStick 0.31
material change diffuse ReviewStick 0.67
material change specular ReviewStick 0.21
proc atlas_overview {pdb base scale ligbonds} {
  set m [mol new $pdb type pdb waitfor all]
  set l [atomselect $m {chain B}]
  if {[$l num] != [llength $ligbonds]} {error "ligand atom count mismatch: $pdb"}
  $l setbonds $ligbonds
  $l delete
  color change rgb blue 0.14 0.39 0.64
  color change rgb silver 0.63 0.76 0.87
  mol ssrecalc $m
  mol delrep 0 $m
  mol representation NewCartoon 0.29 18 4.1 0
  mol selection {protein and not chain B}
  mol material ReviewCartoon
  mol color ColorID 0
  mol addrep $m
  mol representation QuickSurf 1.0 0.65 0.6 1.7
  mol selection {protein and not chain B}
  mol material ReviewEnvelope
  mol color ColorID 6
  mol addrep $m
  mol representation Licorice 0.30 16.0 16.0
  mol selection {chain B}
  mol material ReviewStick
  mol color Element
  mol addrep $m
  display resetview
  set zoom [list [list $scale 0 0 0] [list 0 $scale 0 0] [list 0 0 $scale 0] {0 0 0 1}]
  molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} [list [transidentity] [transidentity] $zoom [transidentity]]
  render Tachyon $base.dat true
  if {![file exists $base.dat] || [file size $base.dat] < 1000} {error "missing Tachyon scene $base"}
  mol delete $m
  puts "ATLAS_RENDERED [file tail $base]"
}
