[hs1]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 1 0'
  length = 0.5
  n_elems = 25

  inner_radius = 0.1
  names = 'region1'
  widths = '0.1'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = 300
[]

[hs2]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 1 0'
  length = '0.5 0.5'
  n_elems = '25 25'
  axial_region_names = 'axregion1 axregion2'

  inner_radius = 0.5
  names = 'region1'
  widths = '0.1'
  n_part_elems = '5'
  solid_properties = 'hs_mat'
  solid_properties_T_ref = '300'

  initial_T = 1000
[]

[hs_coupler]
  type = HeatStructure2DRadiationCouplerRZ
  primary_heat_structure = hs1
  secondary_heat_structure = hs2
  primary_boundary = hs1:outer
  secondary_boundary = hs2:axregion1:inner
  primary_emissivity = ${emissivity1}
  secondary_emissivity = ${emissivity2}
[]
