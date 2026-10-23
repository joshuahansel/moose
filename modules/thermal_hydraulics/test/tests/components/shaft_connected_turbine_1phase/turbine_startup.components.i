[ch_in]
  type = FlowChannel1Phase
  position = '-1 0 0'
  orientation = '1 0 0'
  length = 1
  A = 0.1
  D_h = 1
  fp = eos
[]
[inlet]
  type = InletMassFlowRateTemperature1Phase
  input = 'ch_in:in'
  m_dot = 0
  T = ${T_in}
[]
[turbine]
  type = ShaftConnectedTurbine1Phase
  inlet = 'ch_in:out'
  outlet = 'ch_out:in'
  position = '0 0 0'
  scaling_factor_rhoEV = 1e-5
  A_ref = 0.1
  volume = 0.0002
  inertia_coeff = '1 1 1 1'
  inertia_const = 1.61397
  speed_cr_I = 1e12
  speed_cr_fr = 0
  tau_fr_coeff = '0 0 0 0'
  tau_fr_const = 0
  omega_rated = ${omega_rated}
  D_wheel = 0.4
  head_coefficient = head
  power_coefficient = power
[]
[ch_out]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  A = 0.1
  D_h = 1
  fp = eos
[]
[outlet]
  type = Outlet1Phase
  input = 'ch_out:out'
  p = ${p_out}
[]

[dyno]
  type = ShaftConnectedMotor
  inertia = 10
  torque = dyno_torque_fn
[]
[shaft]
  type = Shaft
  connected_components = 'turbine dyno'
  initial_speed = ${omega_rated}
[]
