# Junction between 3 pipes, 1 of which goes to a dead-end. In the steady-state,
# no flow should go into the dead-end pipe.

[GlobalParams]
  gravity_vector = '0 0 0'

  scaling_factor_1phase = '1 1 1e-5'

  initial_T = 250
  initial_p = 1e5
  initial_vel_x = 1
  initial_vel_y = 0
  initial_vel_z = 0

  closures = simple_closures
[]

[AuxVariables]
  [p0]
    block = 'inlet_pipe outlet_pipe deadend_pipe'
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [p0_kernel]
    type = StagnationPressureAux
    variable = p0
    fp = eos
    e = e
    v = v
    vel = vel
  []
[]

[FluidProperties]
  [eos]
    type = StiffenedGasFluidProperties
    gamma = 1.4
    cv = 725
    q = 0
    q_prime = 0
    p_inf = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [T0]
    type = ParsedFunction
    expression = 'if (x < 1, 300 + 50 * sin(2*pi*x + 1.5*pi), 250)'
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  solve_type = 'NEWTON'
  line_search = 'basic'
  nl_rel_tol = 0
  nl_abs_tol = 1e-6
  nl_max_its = 10

  l_tol = 1e-6
  l_max_its = 10

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  start_time = 0
  end_time = 5
  dt = 0.1

  abort_on_solve_fail = true
[]

[Postprocessors]
  # These post-processors are used for testing that the stagnation pressure in
  # the dead-end pipe is equal to the inlet stagnation pressure.
  [p0_inlet]
    type = SideAverageValue
    variable = p0
    boundary = inlet_pipe:in
  []
  [p0_deadend]
    type = SideAverageValue
    variable = p0
    boundary = deadend_pipe:out
  []
  [test_rel_err]
    type = RelativeDifferencePostprocessor
    value1 = p0_deadend
    value2 = p0_inlet
  []
[]

[Outputs]
  [out]
    type = CSV
    show = test_rel_err
    sync_only = true
    sync_times = '1 2 3 4 5'
  []
  velocity_as_vector = false
[]
