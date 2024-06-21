!include part_hp_base.i

n_radial_vapor = 1
n_azimuthal_hp = 10

[Mesh]
  [circle]
    type = ConcentricCircleMeshGenerator
    has_outer_square = false
    preserve_volumes = true
    rings = '${n_radial_vapor} ${n_radial_clad}'
    num_sectors = ${n_azimuthal_hp}
    radii = '${R_vapor} ${R_clad}'
  []
  [hp]
    type = MeshExtruderGenerator
    input = circle
    extrusion_vector = '0 0 ${fparse length_evap + length_cond}'
    num_layers = ${fparse n_elems_evap + n_elems_cond}
    bottom_sideset = hp_bottom
    top_sideset = hp_top
  []
  [rename_outer]
    type = RenameBoundaryGenerator
    input = hp
    old_boundary = outer
    new_boundary = hp_outer
  []
  [condenser_boundary]
    type = ParsedGenerateSideset
    input = rename_outer
    combinatorial_geometry = 'z > 0.05'
    included_boundaries = 'hp_outer'
    new_sideset_name = condenser
  []
  [hp_coupled_boundary]
    type = ParsedGenerateSideset
    input = condenser_boundary
    combinatorial_geometry = 'z < 0.05'
    included_boundaries = 'hp_outer'
    new_sideset_name = hp_coupled_boundary
  []
  [rename_hp_block]
    type = RenameBlockGenerator
    input = hp_coupled_boundary
    old_block = '1 2'
    new_block = 'vapor clad'
  []
[]

[Materials]
  [vapor_mat]
    type = ADGenericConstantMaterial
    block = vapor
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho_vapor} ${cp_vapor} ${k_vapor}'
  []
  [clad_mat]
    type = ADGenericConstantMaterial
    block = clad
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho_clad} ${cp_clad} ${k_clad}'
  []
[]

[FunctorMaterials]
  [condenser_heat_flux_fmat]
    type = ADConvectionHeatFluxFunctorMaterial
    heat_flux_name = condenser_heat_flux
    htc = ${htc_hx}
    T_fluid = ${T_hx}
    T_solid = temperature
  []
[]

[BCs]
  [condenser_bc]
    type = FunctorNeumannBC
    variable = temperature
    boundary = condenser
    functor = condenser_heat_flux
    flux_is_inward = false
  []
[]

[Postprocessors]
  [heat_rate_hx]
    type = ADSideIntegralFunctorPostprocessor
    boundary = condenser
    functor = condenser_heat_flux
    functor_argument = qp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
