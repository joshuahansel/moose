n_radial_mono = 2
n_azimuthal_mono = 40
n_axial_mono = 5
R_mono_inner = 0.011
R_mono_outer = 0.015
length_mono = 0.05
mono_id = 100

rho_mono = 8000
cp_mono = 500
k_mono = 15

power = 200
S_mono_outer = ${fparse 2 * pi * R_mono_outer * length_mono}
heat_flux_avg = ${fparse power / S_mono_outer}

[Mesh]
  [ring]
    type = AnnularMeshGenerator
    rmin = ${R_mono_inner}
    rmax = ${R_mono_outer}
    nt = ${n_azimuthal_mono}
    nr = ${n_radial_mono}
    quad_subdomain_id = ${mono_id}
  []
  [mono]
    type = MeshExtruderGenerator
    input = ring
    extrusion_vector = '0 0 ${length_mono}'
    num_layers = ${n_axial_mono}
    bottom_sideset = mono_bottom
    top_sideset = mono_top
  []
  [rename_mono_boundary]
    type = RenameBoundaryGenerator
    input = mono
    old_boundary = 'rmin rmax'
    new_boundary = 'mono_inner mono_outer'
  []
  [rename_mono_block]
    type = RenameBlockGenerator
    input = rename_mono_boundary
    old_block = ${mono_id}
    new_block = mono
  []
  [rename_mono_boundary_id]
    type = RenameBoundaryGenerator
    input = rename_mono_block
    old_boundary = '0 1 4 5'
    new_boundary = '200 201 202 203'
  []
[]

[Materials]
  [mono_mat]
    type = ADGenericConstantMaterial
    block = mono
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho_mono} ${cp_mono} ${k_mono}'
  []
[]

[Functions]
  [heat_flux_fn]
    type = ParsedFunction
    expression = '${heat_flux_avg} * (1 + 0.2*cos(3*atan2(y,x)))'
  []
[]

[BCs]
  [heat_bc]
    type = FunctorNeumannBC
    variable = temperature
    boundary = mono_outer
    functor = heat_flux_fn
  []
[]

[Postprocessors]
  [heat_rate_heat]
    type = ADSideIntegralFunctorPostprocessor
    boundary = mono_outer
    functor = heat_flux_fn
    functor_argument = qp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
