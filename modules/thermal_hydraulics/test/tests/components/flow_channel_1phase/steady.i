# Tests that a flow channel can run with Steady executioner.
#
# Note that this solve may fail to converge based on initial guess. For example,
# having a guess with velocity set to zero will fail to converge.

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 2
    T = 500
  []

  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 1 0'
    gravity_vector = '0 0 0'
    length = 1.0
    n_elems = 50
    A = 1.0

    initial_T = 300
    initial_p = 1e5
    initial_vel = 1

    f = 10.0
    closures = simple_closures
    fp = fp

    scaling_factor_1phase = '1 1 1e-5'
  []

  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 2e5
  []

#   [pipe3d]
#     type = VisualizationComponent1Dto3D
#     position = '0 0 0'
#     orientation = '1 0 0'
#     length = 1.0
#     n_elems = 50
#     n_azimuthal_elems = 16
#     n_radial_elems = 5
#     radial_region_widths = 0.2
#     radial_region_names = blah
#   []
[]

# [FunctorMaterials]
#   [fmat]
#     type = ParsedFunctorMaterial
#     expression = '300'
#     # functor_symbols = 'temp'
#     # functor_names = 'T'
#     property_name = booboo
#   []
# []

# [AuxVariables]
#   [yblah]
#   []
# []

# [AuxKernels]
#   [yblah_kernels]
#     type = FunctorAux
#     variable = yblah
#     block = pipe3d:blah
#     functor = booboo
#     execute_on = 'TIMESTEP_END'
#   []
# []

# [Problem]
#   kernel_coverage_check = SKIP_LIST
#   kernel_coverage_block_list = 'pipe3d:blah'
#   material_coverage_check = SKIP_LIST
#   material_coverage_block_list = 'pipe3d:blah'
# []

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady

  solve_type = NEWTON
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7
  nl_max_its = 15

  l_tol = 1e-3
  l_max_its = 10

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  exodus = true
[]
