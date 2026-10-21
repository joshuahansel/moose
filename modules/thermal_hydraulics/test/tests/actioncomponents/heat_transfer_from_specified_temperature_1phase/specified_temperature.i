# Tests a specified wall temperature applied to a flow channel ActionComponent through a native
# HeatTransferFromSpecifiedTemperature1Phase, using the same ClosuresBase-derived object the flow
# channel uses to supply the wall heat transfer coefficient.

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
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 2
    T = 500
  []

  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1.0
    n_elems = 50

    A = 1.0
    fp = fp

    initial_T = 300
    initial_p = 1e5
    initial_vel = 1

    f = 1.0
    closures = simple_closures

    scaling_factor_1phase = '1 1 1e-5'
  []

  [ht_pipe]
    type = HeatTransferFromSpecifiedTemperature1Phase
    flow_channel = pipe
    T_wall = 550
    Hw = 1.0e3
    P_hf = 3.14
  []

  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 2e5
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady

  solve_type = NEWTON
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  exodus = true
[]
