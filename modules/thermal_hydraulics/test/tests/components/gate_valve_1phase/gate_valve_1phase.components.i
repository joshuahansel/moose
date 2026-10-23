[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe1:in'
  m_dot = ${m_dot}
  T = ${T}
[]

[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1.0
  n_elems = 50
  A = ${A}
[]

[volume_junction]
  type = VolumeJunction1Phase
  position = '1 0 0'
  connections = 'pipe1:out pipe2A:in pipe3A:in'
  volume = ${V_junction}
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0
[]

[pipe2A]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '0 1 0'
  length = 0.5
  n_elems = 25
  A = ${A}
[]

[pipe2_valve]
  type = GateValve1Phase
  connections = 'pipe2A:out pipe2B:in'
  open_area_fraction = 0 # (controlled via 'pipe2_valve_control')
[]

[pipe2B]
  type = FlowChannel1Phase
  position = '1 0.5 0'
  orientation = '0 1 0'
  length = 0.5
  n_elems = 25
  A = ${A}
[]

[pipe2_outlet]
  type = Outlet1Phase
  input = 'pipe2B:out'
  p = ${p}
[]

[pipe3A]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 25
  A = ${A}
[]

[pipe3_valve]
  type = GateValve1Phase
  connections = 'pipe3A:out pipe3B:in'
  open_area_fraction = 0 # (controlled via 'pipe3_valve_control')
[]

[pipe3B]
  type = FlowChannel1Phase
  position = '1.5 0 0'
  orientation = '1 0 0'
  length = 0.5
  n_elems = 25
  A = ${A}
[]

[pipe3_outlet]
  type = Outlet1Phase
  input = 'pipe3B:out'
  p = ${p}
[]
