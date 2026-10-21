  [hs]
    type = HeatStructurePlate
    position = '0 0 0'
    orientation = '1 0 0'

    length = 1
    n_elems = 5
    names = 'blk'
    widths = '1'
    depth = 0.5
    n_part_elems = '5'

    initial_T = 350
  []

  [start]
    type = HSBoundarySpecifiedTemperature
    hs = hs
    boundary = hs:start
    T = 300
  []

  [end]
    type = HSBoundarySpecifiedTemperature
    hs = hs
    boundary = hs:end
    T = 400
  []
