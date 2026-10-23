[compressor]
  type = ShaftConnectedCompressor1Phase
  inlet = 'pipe:out'
  outlet = 'pipe:in'
  position = '0 0 0'
  scaling_factor_rhoEV = 1e-5
  volume = ${fparse area*0.45}
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  speed_cr_I = 1e12
  speed_cr_fr = 0
  tau_fr_coeff = '0 0 0 0'
  tau_fr_const = 0
  omega_rated = 200
  mdot_rated = 21.74
  rho0_rated = 1.1812
  c0_rated = 340
  speeds = '0.0 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 2'
  Rp_functions = 'Rp00 Rp04 Rp05 Rp06 Rp07 Rp08 Rp09 Rp10 Rp11 Rp11'
  eff_functions = 'eff00 eff04 eff05 eff06 eff07 eff08 eff09 eff10 eff11 eff11'
[]
[pipe]
  type = FlowChannel1Phase
  position = '0.1 0 0'
  orientation = '1 0 0'
  length = 10
  n_elems = 20
[]

[motor]
  type = ShaftConnectedMotor
  inertia = 1e2
  torque = 100
[]

[shaft]
  type = Shaft
  connected_components = 'motor compressor'
  initial_speed = 100
[]
