[pipe]
  type = ElbowPipe1Phase
  # geometry
  position = '0 0 0'
  orientation = '1 0 0'
  start_angle = 270
  end_angle = 360
  radius = 0.25
  n_elems = 50

  # d = 0.1 m
  A   = 7.8539816340e-03
  D_h  = 1.0000000000e-01
  f = 0.1

  fp = fp
[]

[inlet]
  type = SolidWall1Phase
  input = 'pipe:in'
[]

[outlet]
  type = SolidWall1Phase
  input = 'pipe:out'
[]
