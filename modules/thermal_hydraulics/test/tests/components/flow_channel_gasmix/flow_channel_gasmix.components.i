[pipe]
  type = FlowChannelGasMix
  position = '0 0 0'
  orientation = '1 0 0'
  length = 10.0
  n_elems = 50
  A = 0.2

  initial_mass_fraction = initial_mass_fraction_fn
  initial_p = ${initial_p}
  initial_T = ${initial_T}
  initial_vel = 0

  fp = fp_mix
  closures = 'closures'

  scaling_factor_rhoEA = 1e-5
[]
[inlet]
  type = SolidWallGasMix
  input = 'pipe:in'
[]
[outlet]
  type = SolidWallGasMix
  input = 'pipe:out'
[]
