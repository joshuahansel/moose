n_radial_cylinder = 2
n_azimuthal_cylinder = 10
n_axial_cylinder = 5
R_cylinder = 0.01
length_cylinder = 0.05

[Mesh]
  [circle]
    type = ConcentricCircleMeshGenerator
    has_outer_square = false
    preserve_volumes = true
    rings = ${n_radial_cylinder}
    num_sectors = ${n_azimuthal_cylinder}
    radii = ${R_cylinder}
  []
  [cylinder]
    type = MeshExtruderGenerator
    input = circle
    extrusion_vector = '0 0 ${length_cylinder}'
    num_layers = ${n_axial_cylinder}
    bottom_sideset = cylinder_bottom
    top_sideset = cylinder_top
  []
  [rename_cylinder_boundary]
    type = RenameBoundaryGenerator
    input = cylinder
    old_boundary = outer
    new_boundary = cylinder_outer
  []
  [rename_cylinder_block]
    type = RenameBlockGenerator
    input = rename_cylinder_boundary
    old_block = 1
    new_block = cylinder
  []
[]

[Functions]
  [cylinder_ic_fn]
    type = ParsedFunction
    expression = 'if(abs(sqrt(x^2 + y^2) - ${R_cylinder}) < 1e-3, 500 + 50*cos(5*atan2(y,x)), 500)'
  []
[]

[ICs]
  [cylinder_ic]
    type = FunctionIC
    variable = temperature
    block = cylinder
    function = cylinder_ic_fn
  []
[]
