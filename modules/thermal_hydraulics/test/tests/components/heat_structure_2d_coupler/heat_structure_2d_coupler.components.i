[hs1]
  type = HeatStructureCylindrical
  position = '-0.5 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 5

  names = 'region1'
  widths = '0.5'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = 500
[]

[hs2]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'
  length = '0.5 0.5'
  n_elems = '5 5'
  axial_region_names = 'axregion1 axregion2'

  names = 'region1 region2'
  widths = '0.5 0.2'
  n_part_elems = '5 3'
  solid_properties = 'hs_mat hs_mat'
  solid_properties_T_ref = '300 300'

  initial_T = 300
[]

[hs3]
  type = HeatStructureCylindrical
  position = '0.5 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 5

  names = 'region1'
  widths = '0.5'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = 500
[]

[hs_coupling_1_2]
  type = HeatStructure2DCoupler
  primary_heat_structure = hs2
  secondary_heat_structure = hs1
  primary_boundary = hs2:region1:start
  secondary_boundary = hs1:end
  heat_transfer_coefficient = 1000
[]

[hs_coupling_2_3]
  type = HeatStructure2DCoupler
  primary_heat_structure = hs2
  secondary_heat_structure = hs3
  primary_boundary = hs2:axregion2:outer
  secondary_boundary = hs3:inner
  heat_transfer_coefficient = 500
[]
