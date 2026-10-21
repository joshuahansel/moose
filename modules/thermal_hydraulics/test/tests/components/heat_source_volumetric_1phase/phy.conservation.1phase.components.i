[flow_channel]
  type = FlowChannel1Phase
  position = '0 0 0'
  orientation = '1 0 0'
  length = 1
  n_elems = 10

  A = 1
  f = 0
  fp = fp
  closures = simple_closures

  initial_T = 310
  initial_p = 1e5
  initial_vel = 0
[]
[wall1]
  type = SolidWall1Phase
  input = flow_channel:in
[]
[wall2]
  type = SolidWall1Phase
  input = flow_channel:out
[]

[heat_source]
  type = HeatSourceVolumetric1Phase
  flow_channel = flow_channel
  q = 1e3
[]
