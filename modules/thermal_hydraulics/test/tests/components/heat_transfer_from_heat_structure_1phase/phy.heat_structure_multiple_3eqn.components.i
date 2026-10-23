[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 10
  A = 1
  f = 0
  fp = fp
[]

[hs1]
  type = HeatStructurePlate
  position = '0 -1 0'
  orientation = '1 0 0'
  length = 1
  depth = 1
  n_elems = 10

  solid_properties = 'hs1_mat'
  solid_properties_T_ref = '300'
  n_part_elems = '5'
  widths = '1'
  names = 'solid'

  initial_T = ${T1}
[]

[hs2]
  type = HeatStructurePlate
  position = '0 -1 0'
  orientation = '1 0 0'
  length = 1
  depth = 1
  n_elems = 10

  solid_properties = 'hs2_mat'
  solid_properties_T_ref = '300'
  n_part_elems = '5'
  widths = '1'
  names = 'solid'

  initial_T = ${T2}
[]

[ht1]
  type = HeatTransferFromHeatStructure1Phase
  hs = hs1
  hs_side = outer
  flow_channel = pipe
  Hw = 1e5
  P_hf = 0.5
[]

[ht2]
  type = HeatTransferFromHeatStructure1Phase
  hs = hs2
  hs_side = outer
  flow_channel = pipe
  Hw = 1e5
  P_hf = 0.5
[]

[left]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[right]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
