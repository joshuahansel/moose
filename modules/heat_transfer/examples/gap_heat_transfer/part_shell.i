n_radial_shell = 2
n_azimuthal_shell = 40
n_axial_shell = 5
R_shell_inner = 0.011
R_shell_outer = 0.015
length_shell = 0.05
shell_id = 100

[Mesh]
  [ring]
    type = AnnularMeshGenerator
    rmin = ${R_shell_inner}
    rmax = ${R_shell_outer}
    nt = ${n_azimuthal_shell}
    nr = ${n_radial_shell}
    quad_subdomain_id = ${shell_id}
  []
  [shell]
    type = MeshExtruderGenerator
    input = ring
    extrusion_vector = '0 0 ${length_shell}'
    num_layers = ${n_axial_shell}
    bottom_sideset = shell_bottom
    top_sideset = shell_top
  []
  [rename_shell_boundary]
    type = RenameBoundaryGenerator
    input = shell
    old_boundary = 'rmin rmax'
    new_boundary = 'shell_inner shell_outer'
  []
  [rename_shell_block]
    type = RenameBlockGenerator
    input = rename_shell_boundary
    old_block = ${shell_id}
    new_block = shell
  []
  [rename_shell_boundary_id]
    type = RenameBoundaryGenerator
    input = rename_shell_block
    old_boundary = '0 1 4 5'
    new_boundary = '200 201 202 203'
  []
[]

[ICs]
  [shell_ic]
    type = ConstantIC
    variable = temperature
    block = shell
    value = 500
  []
[]
