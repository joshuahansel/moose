[inlet_bc]
  type = InletVelocityTemperature1Phase
  input = 'inlet:in'
  vel = ${vel}
  T = ${T_in}
[]
[inlet]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 11'
  orientation = '0 0 -1'
  length = 1
  A = 5
[]
[inlet_plenum]
  type = VolumeJunction1Phase
  position = '0 0 10'
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = ${vel}
  K = 100
  connections = 'inlet:out channel1:in channel2:in'
  volume = 1
[]
[channel1]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 10'
  orientation = '0 0 -1'
  length = 10
  A = 4
  D_h = 1
[]
[channel2]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 10'
  orientation = '0 0 -1'
  length = 10
  A = 1
  D_h = 1
[]
[outlet_plenum]
  type = VolumeJunction1Phase
  position = '0 0 0'
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = ${vel}
  connections = 'channel1:out channel2:out outlet:in'
  volume = 1
[]
[outlet]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 0'
  orientation = '0 0 -1'
  length = 1
  A = 5
[]
[outlet_bc]
  type = Outlet1Phase
  p = ${p_out}
  input = 'outlet:out'
[]
