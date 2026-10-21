[SolidProperties]
  [hs_mat]
    type = ThermalFunctionSolidProperties
    rho = 1
    cp = 2
    k = 3
  []
[]

[Executioner]
  type = Transient

  dt = 0.1
  num_steps = 1
[]
