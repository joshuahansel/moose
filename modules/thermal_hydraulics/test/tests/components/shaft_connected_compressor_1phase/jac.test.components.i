[sw1]
  type = SolidWall1Phase
  input = fch1:in
[]

[fch1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
[]

[compressor]
  type = ShaftConnectedCompressor1Phase
  inlet = 'fch1:out'
  outlet = 'fch2:in'
  position = '1 0 0'
  volume = 0.3
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  speed_cr_I = 1e12
  speed_cr_fr = 0
  tau_fr_coeff = '0 0 9.084 0'
  tau_fr_const = 0
  omega_rated = 1476.6
  mdot_rated = 2
  rho0_rated = 1.3
  c0_rated = 350
  speeds = '-1.0 0.0 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.5'
  Rp_functions = 'Rp00 Rp00 Rp04 Rp05 Rp06 Rp07 Rp08 Rp09 Rp10 Rp11 Rp11'
  eff_functions = 'eff00 eff00 eff04 eff05 eff06 eff07 eff08 eff09 eff10 eff11 eff11'
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
  connected_components = 'compressor'
  initial_speed = 1476.6
[]
