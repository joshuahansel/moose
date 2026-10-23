[sw1]
  type = SolidWall1Phase
  input = fch1:in
[]

[fch1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  initial_p = 2e6
[]

[turbine]
  type = ShaftConnectedTurbine1Phase
  inlet = 'fch1:out'
  outlet = 'fch2:in'
  position = '1 0 0'
  volume = 0.3
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  speed_cr_I = 1e12
  speed_cr_fr = 0
  tau_fr_coeff = '0 0 12 0'
  tau_fr_const = 0
  omega_rated = 295
  D_wheel = 0.4
  head_coefficient = head
  power_coefficient = power
[]

[fch2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
[]

[sw2]
  type = SolidWall1Phase
  input = fch2:out
[]

[shaft]
  type = Shaft
  connected_components = 'turbine'
  initial_speed = 300
[]
