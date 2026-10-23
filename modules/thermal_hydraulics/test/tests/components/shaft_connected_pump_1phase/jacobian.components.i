[fch1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
  A = 1
[]

[pump]
  type = ShaftConnectedPump1Phase
  inlet = 'fch1:out'
  outlet = 'fch2:in'
  position = '1 0 0'
  volume = 0.3
  A_ref = 1
  initial_vel_x = 0
  initial_vel_y = 0
  initial_vel_z = 0
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  omega_rated = 314
  speed_cr_I = 1e12
  speed_cr_fr = 0
  torque_rated = 47.1825
  volumetric_rated = 1
  head_rated = 58.52
  tau_fr_coeff = '0 0 9.084 0'
  tau_fr_const = 0
  head = head_fcn
  torque_hydraulic = torque_fcn
  density_rated = 1
[]

[fch2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 2
  A = 1
[]

[shaft]
  type = Shaft
  connected_components = 'pump'
  initial_speed = 1
[]
