# This input file is used to test the gate valve component.

# This problem consists of a T junction of 3 pipes. The inlet pipe is one of the
# 2 pipes of the "top" of the T. The other 2 pipes each have a gate valve.
# Initially, one of the 2 outlet pipes has an open valve and the other has a
# closed valve. Later in the transient, the valves gradually open/close to switch
# the outlet flow direction.

p = 1.0e5
T = 300.0
rho = 1.161430436 # @ 1e5 Pa, 300 K

D = 0.1
A = ${fparse pi * D^2 / 4.0}
V_junction = ${fparse pi * D^3 / 4.0}

vel_in = 2.0
m_dot = ${fparse rho * vel_in * A}

t_begin = 0.3
delta_t_open = 0.1

[GlobalParams]
  gravity_vector = '0 0 0'

  closures = simple_closures
  fp = fp

  f = 0.0

  initial_T = ${T}
  initial_p = ${p}
  initial_vel = 0
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
    molar_mass = 0.02897
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [pipe3_open_fn]
    type = TimeRampFunction
    initial_value = 1
    final_value = 0
    initial_time = ${t_begin}
    ramp_duration = ${delta_t_open}
  []
  [pipe2_open_fn]
    type = ParsedFunction
    expression = '1 - pipe3_phi'
    symbol_names = 'pipe3_phi'
    symbol_values = 'pipe3_open_fn'
  []
[]

[ChainControls]
  [pipe2_valve_control_fn]
    type = GetFunctionValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    function = pipe2_open_fn
  []
  [pipe2_valve_control]
    type = SetRealValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    parameter = component::pipe2_valve/open_area_fraction
    value = pipe2_valve_control_fn:value
  []
  [pipe3_valve_control_fn]
    type = GetFunctionValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    function = pipe3_open_fn
  []
  [pipe3_valve_control]
    type = SetRealValueChainControl
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    parameter = component::pipe3_valve/open_area_fraction
    value = pipe3_valve_control_fn:value
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

  solve_type = PJFNK
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 20

  l_tol = 1e-4

  start_time = 0.0
  end_time = 1.0
  dt = 0.01
  abort_on_solve_fail = true
[]

[Outputs]
  print_linear_residuals = false
  velocity_as_vector = false
  [exodus]
    type = Exodus
    show = 'p T vel'
    output_dimension = 3 # classic components write 3D meshes
    append_object_name = false
  []
  [console]
    type = Console
    max_rows = 1
  []
[]
