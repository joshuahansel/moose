!include part_base.i
!include part_hp_base.i
!include part_htc_gap.i

[AuxVariables]
  [T_mono]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${T_initial}
  []
  [htc_gap]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${htc_gap}
  []
[]

[Components]
  [hp]
    type = HeatStructureCylindrical
    position = '0 0 0'
    orientation = '0 0 1'
    length = '${length_evap} ${length_cond}'
    n_elems = '${n_elems_evap} ${n_elems_cond}'
    axial_region_names = 'evap cond'

    widths = '${R_vapor} ${fparse R_clad - R_vapor}'
    n_part_elems = '1 ${n_radial_clad}'
    names = 'vapor clad'

    initial_T = ${T_initial}
  []
  [hp_evap_bc]
    type = HSBoundaryExternalAppConvection
    hs = hp
    boundary = 'hp:evap:outer'
    T_ext = T_mono
    htc_ext = htc_gap
    add_T_ext = false
    add_htc_ext = false
  []
  [hp_cond_bc]
    type = HSBoundaryAmbientConvection
    hs = hp
    boundary = 'hp:cond:outer'
    T_ambient = ${T_hx}
    htc_ambient = ${htc_hx}
  []
[]

[Materials]
  [vapor_mat]
    type = ADGenericConstantMaterial
    block = hp:vapor
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho_vapor} ${cp_vapor} ${k_vapor}'
  []
  [clad_mat]
    type = ADGenericConstantMaterial
    block = hp:clad
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho_clad} ${cp_clad} ${k_clad}'
  []
[]

[UserObjects]
  [T_hp_uo]
    type = LayeredSideAverage
    variable = T_solid
    boundary = 'hp:evap:outer'
    num_layers = ${n_elems_evap}
    direction = z
    direction_min = 0
    direction_max = ${length_evap}
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [gap_heat_rate_main]
    type = Receiver
  []
[]
