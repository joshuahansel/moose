# heat_rate = 1000

inlet_p0 = ${units 1 atm -> Pa}
inlet_T0 = 300

source_p = 1e5
source_T = 500

# Offsetting from axis by 1e-7 just to test that it still works for some small, nonzero distance
junction_position = '7.1 1e-7 0'

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

[Components]
  [inlet]
    type = InletStagnationPressureTemperature1Phase
    input = 'pipe:in'
    p0 = ${inlet_p0}
    T0 = ${inlet_T0}
  []

  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    gravity_vector = '0 0 0'
    length = 10
    n_elems = 10
    A = 1.0
    fp = fp
    closures = simple_closures
    f = 0

    initial_p = ${inlet_p0}
    initial_T = ${inlet_T0}
    initial_vel = 0

    scaling_factor_1phase = '1 1 1e-5'
  []

  [right_wall]
    type = SolidWall1Phase
    input = 'pipe:out'
  []
[]

[DiracKernels]
  [junction_rhoA]
    type = PointJunction1PhaseDiracKernel
    block = pipe
    variable = rhoA
    pressure = ${source_p}
    temperature = ${source_T}
    point = ${junction_position}
  []
  [junction_rhoEA]
    type = PointJunction1PhaseDiracKernel
    block = pipe
    variable = rhoEA
    pressure = ${source_p}
    temperature = ${source_T}
    point = ${junction_position}
  []
[]

[Postprocessors]
  [total_energy]
    type = ElementIntegralVariablePostprocessor
    variable = rhoEA
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = total_energy
    execute_on = 'INITIAL TIMESTEP_END'
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
  dt = 1.0
  num_steps = 10
  abort_on_solve_fail = true

  solve_type = 'NEWTON'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  nl_max_its = 30

  l_tol = 1e-3
  l_max_its = 100
[]

[Outputs]
  exodus = true
[]
