[pipe1]
  type = FlowChannel1Phase
  position = '0.2 0 0'
  orientation = '0 1 0'
  length = 1
  n_elems = 50

  A   = 9.6858407346e-01
  D_h  = 6.1661977237e+00

  f = 0.01

  fp = eos
[]

[hs]
  type = HeatStructureCylindrical
  position = '0.1 1 0'
  orientation = '0 -1 0'
  length = 1
  n_elems = 50

  solid_properties = 'wall-mat'
  solid_properties_T_ref = '300'
  n_part_elems = 3
  widths = '0.1'
  names = 'wall'

  initial_T = T_init
[]

[hxconn]
  type = HeatTransferFromHeatStructure1Phase
  hs = hs
  hs_side = outer
  flow_channel = pipe1
  Hw = 0
  P_hf = 6.2831853072e-01
[]

[inlet]
  type = SolidWall1Phase
  input = 'pipe1:in'
[]

[outlet]
  type = SolidWall1Phase
  input = 'pipe1:out'
[]
