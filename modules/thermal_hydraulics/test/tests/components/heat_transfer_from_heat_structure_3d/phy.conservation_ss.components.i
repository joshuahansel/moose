[in1]
  type = InletVelocityTemperature1Phase
  input = 'fch1:in'
  vel = 1
  T = 300
[]
[fch1]
  type = FlowChannel1Phase
  position = '0.15 0 0'
  orientation = '0 0 1'
  fp = fp
  n_elems = 10
  length = 1
  initial_T = 300
  initial_p = 1.01e5
  initial_vel = 1
  closures = simple_closures
  A = 0.00314159
  f = 0.0
[]
[out1]
  type = Outlet1Phase
  input = 'fch1:out'
  p = 1.01e5
[]

[in2]
  type = InletVelocityTemperature1Phase
  input = 'fch2:in'
  vel = 1
  T = 350
[]
[fch2]
  type = FlowChannel1Phase
  position = '0 0.15 0'
  orientation = '0 0 1'
  fp = fp
  n_elems = 10
  length = 1
  initial_T = 350
  initial_p = 1.01e5
  initial_vel = 1
  closures = simple_closures
  A = 0.00314159
  f = 0
[]
[out2]
  type = Outlet1Phase
  input = 'fch2:out'
  p = 1.01e5
[]

[blk]
  type = HeatStructureFromFile3D
  file = mesh.e
  position = '0 0 0'
  initial_T = 325
[]
[ht]
  type = HeatTransferFromHeatStructure3D1Phase
  flow_channels = 'fch1 fch2'
  hs = blk
  boundary = blk:rmin
  Hw = 10000
  P_hf = ${P_hf}
[]
