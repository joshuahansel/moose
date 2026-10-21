[hs1]
  type = HeatStructureCylindrical
  position = ${position1}
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'body'
  widths = '${radius_123}'
  n_part_elems = '${n_radial_elems_123}'
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
  widths = '${radius_123}'
  n_part_elems = '${n_radial_elems_123}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T2}
[]
[hs3]
  type = HeatStructureCylindrical
  position = ${position3}
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'body'
  widths = '${radius_123}'
  n_part_elems = '${n_radial_elems_123}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T3}
[]
[hs4]
  type = HeatStructureCylindrical
  position = ${position4}
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  inner_radius = ${inner_radius_4}
  names = 'body'
  widths = '${thickness_4}'
  n_part_elems = '${n_radial_elems_4}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T4}
[]
[hs_coupler]
  type = HSCoupler2D2DRadiation
  heat_structures = 'hs1 hs2 hs3 hs4'
  boundaries = 'hs1:outer hs2:outer hs3:outer hs4:inner'
  emissivities = '${emissivity1} ${emissivity2} ${emissivity3} ${emissivity4}'
  include_environment = false
  view_factors = '
    0 ${view_factor_12} ${view_factor_12} ${view_factor_14};
    ${view_factor_12} 0 ${view_factor_12} ${view_factor_14};
    ${view_factor_12} ${view_factor_12} 0 ${view_factor_14};
    ${view_factor_41} ${view_factor_41} ${view_factor_41} ${view_factor_44}'
[]
