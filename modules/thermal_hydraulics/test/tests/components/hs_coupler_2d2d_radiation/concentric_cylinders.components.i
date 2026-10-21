# Setup with HSCoupler2D2DRadiation

[hs1A]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'region1'
  widths = '${outer_radius1}'
  n_part_elems = '${n_radial_elems1}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T1}
[]
[hs2A]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  inner_radius = ${inner_radius2}
  names = 'region1'
  widths = '${thickness2}'
  n_part_elems = '${n_radial_elems2}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T2}
[]
[hs_couplerA]
  type = HSCoupler2D2DRadiation
  heat_structures = 'hs1A hs2A'
  boundaries = 'hs1A:outer hs2A:inner'
  emissivities = '${emissivity1} ${emissivity2}'
  include_environment = false
  view_factors = '0.0 1.0; ${view_factor21} ${view_factor22}'
[]

# Setup with HeatStructure2DRadiationCouplerRZ

[hs1B]
  type = HeatStructureCylindrical
  position = '0 ${y_shiftB} 0'
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  names = 'region1'
  widths = '${outer_radius1}'
  n_part_elems = '${n_radial_elems1}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T1}
[]
[hs2B]
  type = HeatStructureCylindrical
  position = '0 ${y_shiftB} 0'
  orientation = ${orientation}
  length = ${length}
  n_elems = ${n_axial_elems}

  inner_radius = ${inner_radius2}
  names = 'region1'
  widths = '${thickness2}'
  n_part_elems = '${n_radial_elems2}'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '${T_ref}'

  initial_T = ${initial_T2}
[]
[hs_couplerB]
  type = HeatStructure2DRadiationCouplerRZ
  primary_heat_structure = hs1B
  secondary_heat_structure = hs2B
  primary_boundary = hs1B:outer
  secondary_boundary = hs2B:inner
  primary_emissivity = ${emissivity1}
  secondary_emissivity = ${emissivity2}
[]
