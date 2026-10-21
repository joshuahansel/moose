  [E_tot]
    type = ADHeatStructureEnergyRZ
    block = 'hs1:region1 hs2:region1'
    axis_dir = '1 1 0'
    axis_point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_tot_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = E_tot
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T1]
    type = SideAverageValue
    variable = T_solid
    boundary = hs1:outer
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T2]
    type = SideAverageValue
    variable = T_solid
    boundary = hs2:axregion1:inner
    execute_on = 'INITIAL TIMESTEP_END'
  []
