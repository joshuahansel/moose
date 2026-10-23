[in]
  type = InletStagnationPressureTemperature1Phase
  input = fch1:in
  p0 = 1.1e5
  T0 = 300
[]

[fch1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 1
  A = 1
[]

[junction]
  type = JunctionParallelChannels1Phase
  connections = 'fch1:out fch2:in'
  position = '1 0 0'
  volume = 0.3
  initial_vel_x = 1
  initial_vel_y = 0
  initial_vel_z = 0
[]

[fch2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 1
  A = 1.5
[]

[out]
  type = Outlet1Phase
  input = fch2:out
  p = 1e5
[]
