[pipe1]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  A = 1
  n_elems = 20
[]

[junction1]
  type = JunctionParallelChannels1Phase
  connections = 'pipe1:out pipe2:in'
  position = '1 0 0'
  volume = 1e-2
  K = 0
[]

[pipe2]
  type = FlowChannel1Phase
  position = '1 0 0'
  orientation = '1 0 0'
  length = 1
  A = 0.5
  n_elems = 20
[]

[junction2]
  type = JunctionParallelChannels1Phase
  connections = 'pipe2:out pipe1:in'
  position = '1 0 0'
  volume = 1e-2
[]
