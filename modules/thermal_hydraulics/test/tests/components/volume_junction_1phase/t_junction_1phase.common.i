# Junction between 3 pipes, 1 of which goes to a dead-end. All ends are walls,
# and 1 of the pipes is pressurized higher than the others.

A_big = 1
A_small = 0.5

[GlobalParams]
  gravity_vector = '0 0 0'

  scaling_factor_1phase = '1 1 1e-5'
  scaling_factor_rhoV  = 1
  scaling_factor_rhouV = 1
  scaling_factor_rhovV = 1
  scaling_factor_rhowV = 1
  scaling_factor_rhoEV = 1e-5

  initial_T = 300
  initial_vel = 0

  n_elems = 20
  length = 1

  f = 0

  fp = fp

  rdg_slope_reconstruction = minmod

  closures = simple_closures
[]

[FluidProperties]
  [fp]
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

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  start_time = 0
  end_time = 5
  dt = 0.05
  num_steps = 5
  abort_on_solve_fail = true

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 10

  l_tol = 1e-3
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Postprocessors]
  # mass conservation
  [mass_pipes]
    type = ElementIntegralVariablePostprocessor
    variable = rhoA
    block = 'pipe1 pipe2 pipe3'
    execute_on = 'initial timestep_end'
  []
  [mass_junction]
    type = ElementAverageValue
    variable = rhoV
    block = 'junction'
    execute_on = 'initial timestep_end'
  []
  [mass_tot]
    type = SumPostprocessor
    values = 'mass_pipes mass_junction'
    execute_on = 'initial timestep_end'
  []
  [mass_tot_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = mass_tot
    compute_relative_change = true
    execute_on = 'initial timestep_end'
  []

  # energy conservation
  [E_pipes]
    type = ElementIntegralVariablePostprocessor
    variable = rhoEA
    block = 'pipe1 pipe2 pipe3'
    execute_on = 'initial timestep_end'
  []
  [E_junction]
    type = ElementAverageValue
    variable = rhoEV
    block = 'junction'
    execute_on = 'initial timestep_end'
  []
  [E_tot]
    type = SumPostprocessor
    values = 'E_pipes E_junction'
    execute_on = 'initial timestep_end'
  []
  [E_tot_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = E_tot
    compute_relative_change = true
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  [out]
    type = CSV
    show = 'mass_tot_change E_tot_change'
  []
[]
