[fch1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '0 1 0'
  length = 1
  n_elems = 2
  A = 1
  closures = simple_closures
  fp = fp
  f = 0.01

  initial_p = 1e5
  initial_T = 300
  initial_vel = 0
[]

[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '1 0 0'

  length = 1
  n_elems = 2
  names = 'blk'
  widths = '0.1'
  n_part_elems = '1'
  solid_properties = 'mat'
  solid_properties_T_ref = '300'

  initial_T = 300
[]

[hx]
  type = HeatTransferFromHeatStructure1Phase
  hs = hs
  hs_side = START
  flow_channel = fch1
  Hw = 0
[]

[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'fch1:in'
  m_dot = 1
  T = 300
[]

[outlet]
  type = Outlet1Phase
  input = 'fch1:out'
  p = 1e5
[]
