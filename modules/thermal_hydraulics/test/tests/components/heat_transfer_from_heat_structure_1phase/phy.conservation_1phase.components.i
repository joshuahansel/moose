[left_wall]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[pipe]
  type = FlowChannel1Phase
  fp = fp

  position = '0 2 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 5
  A = 1.0

  initial_T = 300
  initial_p = 1e5
  initial_vel = 0

  f = 0
[]

[right_wall]
  type = SolidWall1Phase
  input = 'pipe:out'
[]

[heat_transfer]
  type = HeatTransferFromHeatStructure1Phase
  flow_channel = pipe
  hs = heat_structure
  hs_side = inner
  Hw = 1e3
[]

[heat_structure]
  #type = set externally
  num_rods = 5

  position = '0 2 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 5

  names = 'main'
  solid_properties = 'main-material'
  solid_properties_T_ref = '300'
  widths = '1.0'
  n_part_elems = '5'

  initial_T = T0_fn
[]
