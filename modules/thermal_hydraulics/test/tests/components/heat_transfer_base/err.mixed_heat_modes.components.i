[pipe]
  type = FlowChannel1Phase
  fp = fp_water
  position = '0 0 0'
  orientation = '1 0 0'
  A = 1
  f = 0
  length = 1
  n_elems = 1
[]

[ht1]
  type = HeatTransferFromHeatFlux1Phase
  flow_channel = pipe
  q_wall = 1
  P_hf = 1
  Hw = 1
[]

[ht2]
  type = HeatTransferFromSpecifiedTemperature1Phase
  flow_channel = pipe
  T_wall = 300
  P_hf = 1
  Hw = 1
[]

[left]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[right]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
