# Verifies that connecting more than one heat transfer component to a single flow channel
# ActionComponent is an error (only one connection is currently supported).

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[ActionComponents]
  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1.0
    n_elems = 10

    A = 1.0
    fp = fp

    initial_T = 300
    initial_p = 1e5
    initial_vel = 1

    f = 1.0
    closures = simple_closures

    scaling_factor_1phase = '1 1 1e-5'
  []

  [ht1]
    type = HeatTransferFromHeatFlux1Phase
    flow_channel = pipe
    q_wall = 1.0e4
    Hw = 1.0e4
    P_hf = 3.14
  []

  [ht2]
    type = HeatTransferFromHeatFlux1Phase
    flow_channel = pipe
    q_wall = 2.0e4
    Hw = 1.0e4
    P_hf = 3.14
  []
[]

[Executioner]
  type = Steady
[]
