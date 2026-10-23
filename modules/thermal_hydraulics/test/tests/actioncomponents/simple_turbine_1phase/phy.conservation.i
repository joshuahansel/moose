!include ../../components/simple_turbine_1phase/phy.conservation.common.i

[ActionComponents]
  !include ../../components/simple_turbine_1phase/phy.conservation.components.i
[]

[Postprocessors]
  !include phy.conservation.vel.i
[]
