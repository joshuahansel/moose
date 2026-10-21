[hs1]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'
  length = '${L_hs}'
  n_elems = '5'

  names = 'region1'
  widths = '${R1}'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = ${initial_T1}
[]

[hs2]
  type = HeatStructureCylindrical
  position = '0 0.3 0'
  orientation = '1 0 0'
  length = '${L_hs}'
  n_elems = '5'

  names = 'region1'
  widths = '${R2}'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = ${initial_T2}
[]

[hs_coupling]
  type = HeatStructure2DCoupler
  primary_heat_structure = hs1
  secondary_heat_structure = hs2
  primary_boundary = hs1:outer
  secondary_boundary = hs2:outer
  heat_transfer_coefficient = ${htc}
[]
