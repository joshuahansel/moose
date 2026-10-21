  [E_hs1]
    type = ADHeatStructureEnergyRZ
    block = 'hs1:body'
    axis_dir = ${orientation}
    axis_point = ${position1}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs2]
    type = ADHeatStructureEnergyRZ
    block = 'hs2:body'
    axis_dir = ${orientation}
    axis_point = ${position2}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs3]
    type = ADHeatStructureEnergyRZ
    block = 'hs3:body'
    axis_dir = ${orientation}
    axis_point = ${position3}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [E_hs4]
    type = ADHeatStructureEnergyRZ
    block = 'hs4:body'
    axis_dir = ${orientation}
    axis_point = ${position4}
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
