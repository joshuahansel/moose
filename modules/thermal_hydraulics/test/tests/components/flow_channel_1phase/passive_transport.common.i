A = 0.1
vel_inlet = 10.0
p_outlet = 1e5
rho_inlet = 1.162633159 # rho @ (1e5 Pa, 300 K)

mdot_inlet = ${fparse rho_inlet * vel_inlet * A}
T_inlet = 300

p_initial = ${p_outlet}
T_initial = ${T_inlet}
vel_initial = ${vel_inlet}

[GlobalParams]
  gravity_vector = '0 0 0'
  scaling_factor_1phase = '1 1 1e-5'
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
  [passiveA_ic]
    type = CosineHumpFunction
    axis = x
    hump_begin_value = 0
    hump_center_value = 0.1
    hump_center_position = 50
    hump_width = 40
  []
  [passiveB_ic]
    type = PiecewiseConstant
    axis = x
    x = '0 50'
    y = '0.1 0'
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
  scheme = bdf2

  start_time = 0
  dt = 0.2
  num_steps = 5

  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  nl_max_its = 30

  l_tol = 1e-3
  l_max_its = 10
[]

[Outputs]
  exodus = true
  show = 'passiveA_times_area passiveB_times_area'
[]
