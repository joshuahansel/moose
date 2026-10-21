  [energy_hs2d]
    type = ADHeatStructureEnergy
    block = 'hs2d:pipe'
    plate_depth = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_hs3d]
    type = ADHeatStructureEnergy3D
    block = 'hs3d:0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy]
    type = SumPostprocessor
    values = 'energy_hs3d energy_hs2d'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [energy_change]
    type = ChangeOverTimePostprocessor
    change_with_respect_to_initial = true
    postprocessor = total_energy
    compute_relative_change = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
