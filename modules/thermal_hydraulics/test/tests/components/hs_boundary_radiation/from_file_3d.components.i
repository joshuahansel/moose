[hs]
  type = HeatStructureFromFile3D
  file = box.e
  position = '0 0 0'
  initial_T = ${T_hs}
[]

[hs_boundary]
  type = HSBoundaryRadiation
  boundary = 'hs:left'
  hs = hs
  T_ambient = ${T_ambient}
  emissivity = ${emissivity}
  view_factor = ${view_factor}
  scale = ${scale}
[]
