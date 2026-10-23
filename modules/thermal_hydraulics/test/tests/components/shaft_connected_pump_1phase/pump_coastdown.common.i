# Pump data used in this test comes from the Semiscale Program, summarized in NUREG/CR-4945

initial_T = 393.15
area = 1e-2
dt = 0.005

[GlobalParams]
  initial_p = 1.4E+07
  initial_T = ${initial_T}
  initial_vel = 0.01
  initial_vel_x = 0.01
  initial_vel_y = 0
  initial_vel_z = 0
  A = ${area}
  A_ref = ${area}
  f = 100
  scaling_factor_1phase = '1 1 1e-3'
  closures = simple_closures
  rdg_slope_reconstruction = minmod
  fp = fp
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [head_fcn]
    type = PiecewiseLinear
    data_file = semiscale_head_data.csv
    format = columns
  []
  [torque_fcn]
    type = PiecewiseLinear
    data_file = semiscale_torque_data.csv
    format = columns
  []
[]

[Postprocessors]
  [vel_avg]
    type = ElementAverageValue
    variable = vel
    block = 'pipe'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [hydraulic_torque]
    type = ElementAverageValue
    variable = hydraulic_torque
    block = 'pump'
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  dt = ${dt}
  num_steps = 40

  solve_type = 'NEWTON'
  line_search = 'basic'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 15

  l_tol = 1e-4
  l_max_its = 10

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  velocity_as_vector = false
  file_base = 'pump_coastdown'
  [csv]
    type = CSV
    show = 'shaft:omega vel_avg'
  []
[]
