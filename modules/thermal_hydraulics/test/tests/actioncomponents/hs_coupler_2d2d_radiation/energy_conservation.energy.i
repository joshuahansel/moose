  [E_hs1]
    type = ADHeatStructureEnergy
    block = 'hs1:body'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs2]
    type = ADHeatStructureEnergy
    block = 'hs2:body'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs3]
    type = ADHeatStructureEnergy
    block = 'hs3:body'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs4]
    type = ADHeatStructureEnergy
    block = 'hs4:body'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_tot]
    type = ParsedPostprocessor
    expression = 'E_hs1 + E_hs2 + E_hs3 + E_hs4'
    pp_names = 'E_hs1 E_hs2 E_hs3 E_hs4'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_tot_err]
    type = ChangeOverTimePostprocessor
    postprocessor = E_tot
    take_absolute_value = true
    change_with_respect_to_initial = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
