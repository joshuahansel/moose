# Test dispersive part of PorousFlowDispersiveFlux kernel by setting diffusion
# coefficients to zero. A pressure gradient is applied over the mesh to give a
# uniform velocity. Gravity is set to zero.
# Mass fraction is set to 1 on the left hand side and 0 on the right hand side.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmax = 10
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [./pp]
  [../]
  [./massfrac0]
  [../]
[]

[AuxVariables]
  [./velocity]
    family = MONOMIAL
    order = FIRST
  [../]
[]

[AuxKernels]
  [./velocity]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity
    component = x
  [../]
[]

[ICs]
  [./pp]
    type = FunctionIC
    variable = pp
    function = pic
  [../]
  [./massfrac0]
    type = ConstantIC
    variable = massfrac0
    value = 0
  [../]
[]

[Functions]
  [./pic]
    type = ParsedFunction
    value = 1.1e5-x*1e3
  [../]
[]

[BCs]
  [./xleft]
    type = PresetBC
    value = 1
    variable = massfrac0
    boundary = left
  [../]
  [./xright]
    type = PresetBC
    value = 0
    variable = massfrac0
    boundary = right
  [../]
  [./pright]
    type = PresetBC
    variable = pp
    boundary = right
    value = 1e5
  [../]
  [./pleft]
    type = PresetBC
    variable = pp
    boundary = left
    value = 1.1e5
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  [../]
  [./adv0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
  [../]
  [./diff0]
    type = PorousFlowDispersiveFlux
    variable = pp
    disp_trans = 0
    disp_long = 0.2
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = massfrac0
  [../]
  [./adv1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = massfrac0
  [../]
  [./diff1]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = massfrac0
    disp_trans = 0
    disp_long = 0.2
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp massfrac0'
    number_fluid_phases = 1
    number_fluid_components = 2
  [../]
[]

[Materials]
  [./nnn]
    type = PorousFlowNodeNumber
    on_initial_only = true
  [../]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./ppss]
    type = PorousFlow1PhaseP
    porepressure = pp
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = massfrac0
  [../]
  [./dens0]
    type = PorousFlowDensityConstBulk
    density_P0 = 1000
    bulk_modulus = 1e9
    phase = 0
  [../]
  [./dens_qp_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_qps = true
  [../]
  [./dens_nodal_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density
    include_old = true
  [../]
  [./poro]
    type = PorousFlowPorosityConst
    porosity = 0.3
  [../]
  [./diff]
    type = PorousFlowDiffusionCoeffConst
    diffusion_coeff = '0 0'
    tortuosity = 0.1
  [../]
  [./relp]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  [../]
  [./relp_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_relative_permeability
  [../]
  [./visc0]
    type = PorousFlowViscosityConst
    viscosity = 0.001
    phase = 0
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_viscosity
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-9 0 0 0 1e-9 0 0 0 1e-9'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1e3
  dtmax = 50
  [./TimeStepper]
    type = IterationAdaptiveDT
    growth_factor = 1.5
    cutback_factor = 0.5
    dt = 1
  [../]
[]

[VectorPostprocessors]
  [./xmass]
    type = NodalValueSampler
    sort_by = id
    variable = massfrac0
  [../]
[]

[Outputs]
  csv = true
  execute_on = 'final'
  print_perf_log = true
[]
