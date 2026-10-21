!include ../../components/hs_boundary_ambient_convection/cylindrical.common.i

[ActionComponents]
  !include ../../components/hs_boundary_ambient_convection/cylindrical.components.i
[]

[Postprocessors]
  !include cylindrical.energy.i
[]
