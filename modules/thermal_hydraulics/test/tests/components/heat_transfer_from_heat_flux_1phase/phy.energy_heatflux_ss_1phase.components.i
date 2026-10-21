[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '0 0 1'
  length = 3.66
  n_elems = 10

  A = 1.907720E-04
  D_h = 1.698566E-02
  f = 0.0

  fp = eos
[]

[ht_pipe]
  type = HeatTransferFromHeatFlux1Phase
  flow_channel = pipe
  q_wall = 1.0e5
  Hw = 1.0e4
  P_hf = 4.4925e-2
[]

[inlet]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[outlet]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
