# Tests that a junction that does not support passive transport reports an error if the flow
# channels it connects have passive transport variables.

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
  []
[]

[ActionComponents]
  [pipe1]
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
  [pipe2]
    type = FlowChannel1Phase
    position = '1 0 0'
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
  [junction]
    type = JunctionOneToOne1Phase
    connections = 'pipe1:out pipe2:in'
  []
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe1:in'
    m_dot = 1
    T = 300
    passives = 0.1
  []
  [outlet]
    type = Outlet1Phase
    input = 'pipe2:out'
    p = 1e5
  []
[]

[Executioner]
  type = Steady
[]
