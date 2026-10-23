[pipe1a]
  type = FlowChannel1Phase
  fp = fp
  position = '0 0 0'
  orientation = '1 0 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 0.5
  n_elems = 2
[]

[pipe1b]
  type = FlowChannel1Phase
  fp = fp
  position = '0.5 0 0'
  orientation = '1 0 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 0.5
  n_elems = 2
[]

[pipe2]
  type = FlowChannel1Phase
  fp = fp
  position = '1 0 0'
  orientation = '0 1 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 1
  n_elems = 3
[]

[pipe3]
  type = FlowChannel1Phase
  fp = fp
  position = '1 1 0'
  orientation = '-1 0 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 1
  n_elems = 3
[]

[pipe4]
  type = FlowChannel1Phase
  fp = fp
  position = '0 1 0'
  orientation = '0 -1 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 1
  n_elems = 3
[]

[pipe5]
  type = FlowChannel1Phase
  fp = fp
  position = '1 1 0'
  orientation = '0 1 0'
  A = 0.785398163e-4    #1.0 cm (0.01 m) in diameter, A = 1/4 * PI * d^2
  D_h = 0.01
  f = 0.01
  length = 0.5
  n_elems = 3
[]

[pump]
  type = Pump1Phase
  connections = 'pipe1a:out pipe1b:in'
  head = 1.0
  position = '0.5 0 0'
  volume = 0.785398163e-3
  A_ref = 0.785398163e-4
[]

[junction1]
  type = VolumeJunction1Phase
  connections = 'pipe1b:out pipe2:in'
  volume = 0.785398163e-3
  position = '1 0 0'
[]

[junction2]
  type = VolumeJunction1Phase
  connections = 'pipe2:out pipe3:in pipe5:in'
  volume = 0.785398163e-3
  position = '1 1 0'
[]

[junction3]
  type = VolumeJunction1Phase
  connections = 'pipe3:out pipe4:in'
  volume = 0.785398163e-3
  position = '0 1 0'
[]

[junction4]
  type = VolumeJunction1Phase
  connections = 'pipe4:out pipe1a:in'
  volume = 0.785398163e-3
  position = '0 0 0'
[]

[outlet]
  type = InletStagnationPressureTemperature1Phase
  input = 'pipe5:out'
  p0 = 1.e5
  T0 = 300
[]
