[pipe]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 0'
  orientation = '1 0 0'
  length = 2
  A = 1
  n_elems = 5

  f = 0
[]
[form_loss]
  type = FormLossFromFunction1Phase
  flow_channel = pipe
  K_prime = 0.25
[]
[inlet]
  type = InletDensityVelocity1Phase
  input = 'pipe:in'
  rho = 996.5563397
  vel = 0.5
[]
[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = ${p_out}
[]
