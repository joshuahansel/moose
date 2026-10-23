[pump]
  type = ShaftConnectedPump1Phase
  inlet = 'pipe:out'
  outlet = 'pipe:in'
  position = '0 0 0'
  scaling_factor_rhoEV = 1e-5
  volume = 0.3
  inertia_coeff = '1 1 1 1'
  inertia_const = 0.5
  omega_rated = 314
  speed_cr_I = 1e12
  speed_cr_fr = 0.001
  torque_rated = 47.1825
  volumetric_rated = 1
  head_rated = 58.52
  tau_fr_coeff = '4 0 80 0'
  tau_fr_const = 8
  head = head_fcn
  torque_hydraulic = torque_fcn
  density_rated = 124.2046
[]

[pipe]
  type = FlowChannel1Phase
  position = '0.6096 0 0'
  orientation = '1 0 0'
  length = 10
  n_elems = 20
[]

[shaft]
  type = Shaft
  connected_components = 'pump'
  initial_speed = 1
[]
