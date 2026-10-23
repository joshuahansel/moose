[pipe]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 0'
  orientation = '1 0 0'
  A = 1
  f = 0

  # length and number of elements should be arbitrary for the test
  length = 10
  n_elems = 1
[]

[ht1]
  type = HeatTransferFromHeatFlux1Phase
  flow_channel = pipe
  q_wall = 10e6
  P_hf = 0.2
  Hw = 1
[]

[ht2]
  type = HeatTransferFromHeatFlux1Phase
  flow_channel = pipe
  q_wall = 20e6
  P_hf = 0.4
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
