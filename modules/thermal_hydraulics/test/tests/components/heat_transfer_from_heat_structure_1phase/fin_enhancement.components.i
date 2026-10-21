# pipe1

[pipe1_inlet]
  type = InletMassFlowRateTemperature1Phase
  m_dot = ${mfr_inlet}
  T = ${T_inlet}
  input = 'pipe1:in'
[]
[pipe1]
  type = FlowChannel1Phase
  gravity_vector = '0 0 0'
  position = '0 0 0'
  orientation = '0 0 1'
  length = ${length}
  n_elems = ${n_elems}
  A = ${area}
  initial_T = ${T_inlet}
  initial_p = ${p_outlet}
  initial_vel = ${vel_inlet}
  fp = fp
  closures = simple_closures
  f = 0
  scaling_factor_1phase = '1 1 1e-5'
[]
[pipe1_outlet]
  type = Outlet1Phase
  p = ${p_outlet}
  input = 'pipe1:out'
[]
[ht1]
  type = HeatTransferFromHeatStructure1Phase
  flow_channel = pipe1
  hs = hs1
  hs_side = inner
  Hw = ${htc}
[]
[hs1]
  type = HeatStructureCylindrical
  position = '0 0 0'
  orientation = '0 0 1'
  length = ${length}
  n_elems = ${n_elems}
  inner_radius = ${fparse 0.5 * diam}
  names = 'main'
  solid_properties = 'sp_ss316'
  solid_properties_T_ref = '300'
  widths = '${t_hs}'
  n_part_elems = '${n_elems_radial}'
  initial_T = ${T_initial_hs}
  scaling_factor_temperature = 1e-5
[]

# pipe 2

[pipe2_inlet]
  type = InletMassFlowRateTemperature1Phase
  m_dot = ${mfr_inlet}
  T = ${T_inlet}
  input = 'pipe2:in'
[]
[pipe2]
  type = FlowChannel1Phase
  gravity_vector = '0 0 0'
  position = '0 0.5 0'
  orientation = '0 0 1'
  length = ${length}
  n_elems = ${n_elems}
  A = ${area}
  initial_T = ${T_inlet}
  initial_p = ${p_outlet}
  initial_vel = ${vel_inlet}
  fp = fp
  closures = simple_closures
  f = 0
  scaling_factor_1phase = '1 1 1e-5'
[]
[pipe2_outlet]
  type = Outlet1Phase
  p = ${p_outlet}
  input = 'pipe2:out'
[]
[ht2]
  type = HeatTransferFromHeatStructure1Phase
  flow_channel = pipe2
  hs = hs2
  hs_side = inner
  Hw = ${htc}
  scale = fin_enhancement
[]
[hs2]
  type = HeatStructureCylindrical
  position = '0 0.5 0'
  orientation = '0 0 1'
  length = ${length}
  n_elems = ${n_elems}
  inner_radius = ${fparse 0.5 * diam}
  names = 'main'
  solid_properties = 'sp_ss316'
  solid_properties_T_ref = '300'
  widths = '${t_hs}'
  n_part_elems = '${n_elems_radial}'
  initial_T = ${T_initial_hs}
  scaling_factor_temperature = 1e-5
[]
