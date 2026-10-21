  [E_hs1]
    type = ADHeatStructureEnergyRZ
    block = 'hs1:region1'
    axis_dir = '1 0 0'
    axis_point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [E_hs2]
    type = ADHeatStructureEnergyRZ
    block = 'hs2:region1'
    axis_dir = '1 0 0'
    axis_point = '0 0.3 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [E_hs1_change]
    type = ChangeOverTimePostprocessor
    postprocessor = E_hs1
    change_with_respect_to_initial = true
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [E_hs2_change]
    type = ChangeOverTimePostprocessor
    postprocessor = E_hs2
    change_with_respect_to_initial = true
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [E_tot]
    type = SumPostprocessor
    values = 'E_hs1 E_hs2'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_tot_change]
    type = ChangeOverTimePostprocessor
    postprocessor = E_tot
    change_with_respect_to_initial = true
    compute_relative_change = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
