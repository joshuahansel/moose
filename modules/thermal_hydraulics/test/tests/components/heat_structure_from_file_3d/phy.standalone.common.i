[Materials]
  [mat]
    type = ADGenericConstantMaterial
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '1 1 1'
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 1
  abort_on_solve_fail = true
[]

[Outputs]
  exodus = true
[]
