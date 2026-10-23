# Tests a junction between 2 flow channels of equal area and orientation. A
# sinusoidal density shape is advected to the right and should not be affected
# by the junction; the solution should be identical to the equivalent
# no-junction solution.

[GlobalParams]
  gravity_vector = '0 0 0'

  initial_p = 1e5
  initial_vel = 1

  A = 25

  f = 0

  fp = fp

  scaling_factor_1phase = '0.04 0.04 0.04e-5'

  closures = simple_closures
[]

[FluidProperties]
  [fp]
    type = StiffenedGasFluidProperties
    gamma = 1.4
    cv = 725
    p_inf = 0
    q = 0
    q_prime = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [T0]
    type = CosineHumpFunction
    axis = x
    hump_center_position = 1
    hump_width = 0.5
    hump_begin_value = 250
    hump_center_value = 300
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
  dt = 0.01
  num_steps = 5
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = ' lu'
  line_search = 'basic'
  nl_rel_tol = 0
  nl_abs_tol = 1e-6
  nl_max_its = 10

  l_tol = 1e-3
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Postprocessors]
  [junction_rho]
    type = ElementAverageValue
    variable = rhoV
    block = 'junction'
    execute_on = 'initial timestep_end'
  []
  [junction_rhou]
    type = ElementAverageValue
    variable = rhouV
    block = 'junction'
    execute_on = 'initial timestep_end'
  []
  [junction_rhoE]
    type = ElementAverageValue
    variable = rhoEV
    block = 'junction'
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  [out]
    type = CSV
    execute_scalars_on = 'none'
    execute_on = 'initial timestep_end'
  []
[]
