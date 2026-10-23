# This input file simulates the Sod shock tube using a junction in the middle
# of the domain. The solution should be exactly equivalent to the problem with
# no junction. This test examines the solutions at the junction connections
# and compares them to gold values generated from a version of this input file
# that has no junction.

[GlobalParams]
  gravity_vector = '0 0 0'

  closures = simple_closures
[]

[Functions]
  [p_ic_fn]
    type = PiecewiseConstant
    axis = x
    direction = right
    x = '0.5 1.0'
    y = '1.0 0.1'
  []

  [T_ic_fn]
    type = PiecewiseConstant
    axis = x
    direction = right
    x = '0.5 1.0'
    y = '1.4 1.12'
  []
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
    molar_mass = 11.64024372
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2

  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  nl_max_its = 60

  l_tol = 1e-4

  start_time = 0.0
  dt = 1e-3
  num_steps = 5
  abort_on_solve_fail = true
[]

[Postprocessors]
  [rhoA_left]
    type = SideAverageValue
    variable = rhoA
    boundary = left_channel:out
    execute_on = 'initial timestep_end'
  []
  [rhouA_left]
    type = SideAverageValue
    variable = rhouA
    boundary = left_channel:out
    execute_on = 'initial timestep_end'
  []
  [rhoEA_left]
    type = SideAverageValue
    variable = rhoEA
    boundary = left_channel:out
    execute_on = 'initial timestep_end'
  []
  [rhoA_right]
    type = SideAverageValue
    variable = rhoA
    boundary = right_channel:in
    execute_on = 'initial timestep_end'
  []
  # rhouA_right is added by tests file
  [rhoEA_right]
    type = SideAverageValue
    variable = rhoEA
    boundary = right_channel:in
    execute_on = 'initial timestep_end'
  []

  # This is present to test that junction sidesets work properly
  [p_avg_junction]
    type = SideAverageValue
    boundary = 'junction'
    variable = p
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  csv = true
  show = 'rhoA_left rhouA_left rhoEA_left rhoA_right rhouA_right rhoEA_right'
  execute_on = 'initial timestep_end'
[]
