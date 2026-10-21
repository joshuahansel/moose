[pipe]
  type = FlowChannel1Phase
  position = '0 0.1 0'
  orientation = '0 0 1'
  length = 4
  n_elems = 2

  A = 8.78882e-5
  D_h = 0.01179
  f = 0.01

  fp = fp
[]

[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '0 0 1'
  length = 4
  n_elems = 2

  names = 'fuel'
  widths = '0.1'
  n_part_elems = '1'
  solid_properties = 'fuel-mat'
  solid_properties_T_ref = '300'

  initial_T = 300
[]

[hx]
  type = HeatTransferFromHeatStructure1Phase
  hs = hs
  hs_side = outer
  flow_channel = pipe
  P_hf = 0.029832559676
[]

[inlet]
  type = InletStagnationPressureTemperature1Phase
  input = 'pipe:in'
  p0 = 1e5
  T0 = 300
[]

[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = 1e5
[]
