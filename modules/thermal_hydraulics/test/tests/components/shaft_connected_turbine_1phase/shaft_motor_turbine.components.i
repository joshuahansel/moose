[turbine]
  type = ShaftConnectedTurbine1Phase
  inlet = 'pipe:out'
  outlet = 'pipe:in'
  position = '0 0 0'
  volume = 0.2
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  speed_cr_I = 1e12
  speed_cr_fr = 0
  tau_fr_coeff = '0 0 0 0'
  tau_fr_const = 0
  omega_rated = 100
  D_wheel = 0.4
  head_coefficient = head
  power_coefficient = power
[]
[pipe]
  type = FlowChannel1Phase
  position = '0.1 0 0'
  orientation = '1 0 0'
  length = 10
  n_elems = 20
  initial_p = 2e6
[]

[dyno]
  type = ShaftConnectedMotor
  inertia = 1e2
  torque = -1e3
[]

[shaft]
  type = Shaft
  connected_components = 'dyno turbine'
  initial_speed = 300
[]
