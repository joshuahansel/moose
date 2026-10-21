# This test tests conservation of energy at steady state for 1-phase flow when a
# heat flux is specified. Conservation is checked by comparing the integral of
# the heat flux against the difference of the boundary fluxes.

[GlobalParams]
  initial_p = 7.0e6
  initial_vel = 0
  initial_T = 513

  gravity_vector = '0.0 0.0 0.0'

  closures = simple_closures
[]

[FluidProperties]
  [eos]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    q = -1167e3
    q_prime = 0
    p_inf = 1.e9
    cv = 1816
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Postprocessors]
  [E]
    type = ElementIntegralVariablePostprocessor
    variable = rhoEA
    execute_on = 'initial timestep_end'
  []

  [E_change]
    type = ChangeOverTimePostprocessor
    postprocessor = E
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  abort_on_solve_fail = true
  dt = 1

  solve_type = 'NEWTON'
  line_search = 'basic'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  nl_max_its = 50

  l_tol = 1e-3
  l_max_its = 60

  start_time = 0
  num_steps = 10
[]

[Outputs]
  [out]
    type = CSV
    show = 'E_change'
  []
  [console]
    type = Console
    show = 'E_change'
  []
[]
