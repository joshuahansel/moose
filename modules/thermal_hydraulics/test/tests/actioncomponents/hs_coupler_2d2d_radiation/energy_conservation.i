!include ../../components/hs_coupler_2d2d_radiation/energy_conservation.common.i

[ActionComponents]
  !include ../../components/hs_coupler_2d2d_radiation/energy_conservation.components.i
[]

[Postprocessors]
  !include energy_conservation.energy.i
[]
