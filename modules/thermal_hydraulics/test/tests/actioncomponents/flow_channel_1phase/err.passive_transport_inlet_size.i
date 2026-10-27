# Tests that an error is reported if the number of passives of an inlet does not match the number
# of passive transport variables of the flow channel.

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
  []
[]

[ActionComponents]
  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 2
    A = 1
    fp = fp
    f = 0
    initial_T = 300
    initial_p = 1e5
    initial_vel = 1
    initial_passives = 0.1
    passives_names = passive
  []
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 1
    T = 300
  []
  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 1e5
  []
[]

[Executioner]
  type = Steady
[]
