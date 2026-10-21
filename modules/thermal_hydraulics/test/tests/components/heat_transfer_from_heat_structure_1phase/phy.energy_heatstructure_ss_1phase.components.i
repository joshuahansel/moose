[reactor]
  type = TotalPower
  power = 1e3
[]

[core:pipe]
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

[core:solid]
  type = HeatStructureCylindrical
  position = '0 -0.0071501 0'
  orientation = '0 0 1'
  length = 3.66
  n_elems = 10

  names =  'FUEL GAP CLAD'
  widths = '6.057900E-03  1.524000E-04  9.398000E-04'
  n_part_elems = '5 1 2'
  solid_properties = 'fuel-mat gap-mat clad-mat'
  solid_properties_T_ref = '300 300 300'

  initial_T = 513
[]

[core:hgen]
  type = HeatSourceFromTotalPower
  hs = core:solid
  regions = 'FUEL'
  power = reactor
  power_fraction = 1
[]

[core:hx]
  type = HeatTransferFromHeatStructure1Phase
  flow_channel = core:pipe
  hs   = core:solid
  hs_side = outer
  Hw = 1.0e4
  P_hf = 4.4925e-2
[]

[inlet]
  type = InletDensityVelocity1Phase
  input = 'core:pipe:in'
  rho = 817.382210128610836
  vel = 2.4
[]

[outlet]
  type = Outlet1Phase
  input = 'core:pipe:out'
  p = 7e6
[]
