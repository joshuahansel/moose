[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 10

  A   = 9.6858407346e-01
  D_h = 6.1661977237e+00
  f = 0.01

  fp = eos
[]

[hxconn]
  type = HeatTransferFromExternalAppHeatFlux1Phase
  flow_channel = pipe1
  Hw = 1e3
[]

[inlet]
  type = SolidWall1Phase
  input = 'pipe1:in'
[]

[outlet]
  type = SolidWall1Phase
  input = 'pipe1:out'
[]
