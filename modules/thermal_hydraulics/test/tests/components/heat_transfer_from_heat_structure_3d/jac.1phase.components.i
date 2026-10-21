[fch]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '0 0 1'
  fp = fp
  n_elems = 6
  length = 1
  initial_T = T_init
  initial_p = 1.01e5
  initial_vel = 0
  closures = simple_closures
  A   = 0.00314159
  D_h  = 0.2
  f = 0.01
[]
[in]
  type = InletVelocityTemperature1Phase
  input = 'fch:in'
  vel = 1
  T = 300
[]
[out]
  type = Outlet1Phase
  input = 'fch:out'
  p = 1.01e5
[]
[blk]
  type = HeatStructureFromFile3D
  file = mesh.e
  position = '0 0 0'
  initial_T = T_init
[]
[ht]
  type = HeatTransferFromHeatStructure3D1Phase
  flow_channels = 'fch'
  hs = blk
  boundary = blk:rmin
  Hw = 10000
  P_hf = 0.1564344650402309
[]
