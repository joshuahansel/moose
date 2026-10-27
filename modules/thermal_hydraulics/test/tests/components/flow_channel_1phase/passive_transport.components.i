[pipe]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 100
  n_elems = 100
  A = ${A}

  initial_T = ${T_initial}
  initial_p = ${p_initial}
  initial_vel = ${vel_initial}
  initial_passives = 'passiveA_ic passiveB_ic'
  passives_names = 'passiveA passiveB'

  fp = fp
  closures = simple_closures
  f = 0
[]
[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'pipe:in'
  m_dot = ${mdot_inlet}
  T = ${T_inlet}
  passives = '0.1 0.05'
[]
[outlet]
  type = Outlet1Phase
  input = 'pipe:out'
  p = ${p_outlet}
[]
