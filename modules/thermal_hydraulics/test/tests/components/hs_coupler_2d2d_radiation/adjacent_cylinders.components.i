[hs1]
  type = HeatStructureCylindrical
  position = ${position1}
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'body'
  widths = '${radius}'
  n_part_elems = '${n_radial_elems}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T1}
[]
[hs2]
  type = HeatStructureCylindrical
  position = ${position2}
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'body'
  widths = '${radius}'
  n_part_elems = '${n_radial_elems}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T2}
[]
[hs_coupler]
  type = HSCoupler2D2DRadiation
  heat_structures = 'hs1 hs2'
  boundaries = 'hs1:outer hs2:outer'
  emissivities = '${emissivity1} ${emissivity2}'
  include_environment = true
  T_environment = ${T3}
  view_factors = '
    0 ${view_factor_12} ${view_factor_13};
    ${view_factor_12} 0 ${view_factor_13};
    0 0 1'
[]
