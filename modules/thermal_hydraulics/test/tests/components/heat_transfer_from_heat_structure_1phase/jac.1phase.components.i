[pipe]
  type = FlowChannel1Phase
  position = '0 0.1 0'
  orientation = '0 0 1'
  length = 2
  n_elems = 1

  A = 8.78882e-5
  D_h = 0.01179
  f = 0.01

  fp = fp
[]

[hs]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '0 0 1'
  length = 2
  n_elems = 1

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
  Hw = 100
  P_hf = 0.029832559676
[]
