!include part_base.i
!include part_shell.i

# emissivity_shell = 0.8
# emissivity_cylinder = 0.8
# gap_conductivity = 1e-5

# gap_distance = 0.001
gap_htc = 3000

initial_T_cylinder = 300

[MultiApps]
  [sub]
    type = TransientMultiApp
    app_type = CombinedApp
    input_files = gap_3d_sub.i
    positions = '0 0 0'
    max_procs_per_app = 1
    output_in_position = true
    output_sub_cycles = false
    sub_cycling = false
    catch_up = true
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[FunctorMaterials]
  [gap_heat_flux_fmat]
    type = ADConvectionHeatFluxFunctorMaterial
    heat_flux_name = gap_heat_flux
    htc = ${gap_htc}
    T_solid = temperature
    T_fluid = T_cylinder
  []
[]

[BCs]
  [gap_heat_flux_bc]
    type = FunctorNeumannBC
    variable = temperature
    boundary = shell_inner
    functor = gap_heat_flux
  []
[]

[AuxVariables]
  [T_cylinder]
    initial_condition = ${initial_T_cylinder}
  []
  # [gap_distance]
  #   initial_condition = ${gap_distance}
  # []
[]

[Transfers]
  [T_from_sub]
    type = MultiAppGeneralFieldNearestLocationTransfer
    source_variable = temperature
    variable = T_cylinder
    from_multi_app = sub
    from_boundaries = cylinder_outer
    to_boundaries = shell_inner
    num_nearest_points = 1
    error_on_miss = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Postprocessors]
  [gap_heat_rate]
    type = ADSideIntegralFunctorPostprocessor
    boundary = shell_inner
    functor = gap_heat_flux
    functor_argument = qp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  fixed_point_max_its = 15
  fixed_point_abs_tol = 1e-8
  fixed_point_rel_tol = 1e-8
[]
