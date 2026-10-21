[AuxVariables]
  [power_density]
    family = MONOMIAL
    order = CONSTANT
    block = 'hs:fuel'
  []
[]

[AuxKernels]
  [mock_power_aux]
    type = ConstantAux
    variable = power_density
    value = 1e9
    block = 'hs:fuel'
  []
[]

[SolidProperties]
  [fuel-mat]
    type = ThermalFunctionSolidProperties
    k = 2.5
    cp = 300.
    rho = 1.032e4
  []
[]

[Executioner]
  type = Transient
  dt = 1.e-2
[]
