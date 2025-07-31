# Heat conduction with fixed temperature on left and convection BC on right:
#
#   d/dx(-k dT/dx) = S'''(T)    (0,1)X(0,1)
#   T = T_inf                    x = 0
#   -k dT/dx = htc (T - T_inf)   x = 1
#
# Source is temperature-dependent and is calculated in the child app:
#   S(T) = B - A * (T - T_inf)^2

k = 15.0
htc = 100.0
T_ambient = 300.0

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
[]

[Variables]
  [T_nodal]
  []
  [T_elem]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [S_parent]
  []
[]

[FunctorMaterials]
  [heat_flux_mat_nodal]
    type = ADParsedFunctorMaterial
    expression = 'htc * (T - T_inf)'
    functor_symbols = 'T T_inf htc'
    functor_names = 'T_nodal ${T_ambient} ${htc}'
    property_name = 'heat_flux_nodal'
  []
  [heat_flux_mat_elem]
    type = ADParsedFunctorMaterial
    expression = 'htc * (T - T_inf)'
    functor_symbols = 'T T_inf htc'
    functor_names = 'T_elem ${T_ambient} ${htc}'
    property_name = 'heat_flux_elem'
  []
[]

[Kernels]
  [T_nodal_diff]
    type = FunctionDiffusion
    variable = T_nodal
    function = ${k}
  []
  [T_nodal_source]
    type = CoupledForce
    variable = T_nodal
    v = S_parent
  []
[]

[FVKernels]
  [T_elem_diff]
    type = FVDiffusion
    variable = T_elem
    coeff = ${k}
  []
  [T_elem_source]
    type = FVCoupledForce
    variable = T_elem
    v = S_parent
  []
[]

[BCs]
  [left_bc_nodal]
    type = DirichletBC
    variable = T_nodal
    boundary = left
    value = ${T_ambient}
  []
  [right_bc_nodal]
    type = FunctorNeumannBC
    variable = T_nodal
    boundary = right
    functor = heat_flux_nodal
    flux_is_inward = false
  []
[]

[FVBCs]
  [left_bc_elem]
    type = FVDirichletBC
    variable = T_elem
    boundary = left
    value = ${T_ambient}
  []
  [right_bc_elem]
    type = FVFunctorNeumannBC
    variable = T_elem
    boundary = right
    functor = heat_flux_elem
    factor = -1
  []
[]

[MultiApps]
  [source_app]
    type = FullSolveMultiApp
    positions = '0 0 0'
    input_files = multiapp_fp_iteration_child.i
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [T_to_child]
    type = MultiAppCopyTransfer
    to_multi_app = source_app
    source_variable = T_nodal
    variable = T_child
    execute_on = 'SAME_AS_MULTIAPP'
  []
  [S_from_child]
    type = MultiAppCopyTransfer
    from_multi_app = source_app
    source_variable = S
    variable = S_parent
    execute_on = 'SAME_AS_MULTIAPP'
  []
[]

[Postprocessors]
  [T_nodal_avg]
    type = AverageNodalVariableValue
    variable = T_nodal
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_elem_avg]
    type = ElementAverageValue
    variable = T_elem
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_nodal_change]
    type = NodalExtremeVariableChange
    variable = T_nodal
    change_over = multiapp_fp_iteration
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_elem_change]
    type = ElementExtremeVariableChange
    variable = T_elem
    change_over = multiapp_fp_iteration
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [fp_it]
    type = NumFixedPointIterations
    get_index_instead_of_count = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Convergence]
  [fp_conv]
    type = IterationCountConvergence
    max_iterations = 5
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  multiapp_fixed_point_convergence = fp_conv
[]

[Outputs]
  [console]
    type = Console
    new_row_detection_columns = all
    execute_postprocessors_on = 'INITIAL MULTIAPP_FIXED_POINT_ITERATION_END'
  []
  [out]
    type = CSV
    new_row_detection_columns = all
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_ITERATION_END'
  []
[]
