[left_wall]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[pipe]
  type = FlowChannel1Phase
  position = '5.0 0 0'
  orientation = '1 0 0'
  length = 5.0
  n_elems = 5
  A = 1.0

  initial_T = 300
  initial_p = 1e5
  initial_vel = 0

  f = 0
  fp = fp
  closures = simple_closures

  scaling_factor_1phase = '1 1 1e-5'
[]

[right_wall]
  type = SolidWall1Phase
  input = 'pipe:out'
[]

[heat_transfer]
  type = HeatTransferFromHeatStructure1Phase
  flow_channel = pipe
  hs = heat_structure
  hs_boundary = heat_structure:region2:inner
  Hw = 1e3
[]

[heat_structure]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'
  length = '5.0 5.0'
  n_elems = '5 5'
  axial_region_names = 'region1 region2'

  names = 'main'
  solid_properties = 'sp_ss316'
  solid_properties_T_ref = '500'
  widths = '1.0'
  n_part_elems = '1'
  inner_radius = 1.0

  initial_T = 500

  scaling_factor_temperature = 1e-8
[]
