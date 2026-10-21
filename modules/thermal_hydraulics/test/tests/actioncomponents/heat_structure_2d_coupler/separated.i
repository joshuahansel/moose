!include ../../components/heat_structure_2d_coupler/separated.common.i

[ActionComponents]
  !include ../../components/heat_structure_2d_coupler/separated.components.i
[]

[Postprocessors]
  !include separated.energy.i
[]
