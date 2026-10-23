[pipe]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 0'
  orientation = '1 0 0'
  length = 2
  A = 1
  n_elems = 10

  f = 0
[]
[form_loss]
  type = FormLossFromExternalApp1Phase
  flow_channel = pipe
[]
[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe:in'
  m_dot = 680
  T = 300
[]
[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = 1e5
[]
