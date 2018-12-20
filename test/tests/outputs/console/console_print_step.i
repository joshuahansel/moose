[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Variables]
  [./u]
    family = SCALAR
    order = FIRST
    initial_condition = 10
  [../]
[]

[ScalarKernels]
  [./func]
    type = ParsedODEKernel
    variable = u
    function = '(u - 5)^3'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  # execute_on = 'timestep_end'
  # color = false
  [./console]
    type = Console
    # fit_mode = 40
  [../]
[]
