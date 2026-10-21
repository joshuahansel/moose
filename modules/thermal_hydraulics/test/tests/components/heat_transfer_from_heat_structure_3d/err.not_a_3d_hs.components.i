[fch]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '0 0 1'
  fp = fp
  n_elems = 6
  length = 1
  initial_T = 300
  initial_p = 1.01e5
  initial_vel = 1
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
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '0 0 1'
  widths = 0.1
  inner_radius = 0.1
  length = 1
  n_elems = 6
  n_part_elems = 1
  initial_T = T_init
  solid_properties = 'mat'
  solid_properties_T_ref = '300'
  names = blk
[]
[ht]
  type = HeatTransferFromHeatStructure3D1Phase
  flow_channels = 'fch'
  hs = blk
  boundary = blk:inner
  Hw = 10000
  P_hf = 0.156434465
[]
