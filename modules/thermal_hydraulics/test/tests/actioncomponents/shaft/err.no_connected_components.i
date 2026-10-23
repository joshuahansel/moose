!include ../../components/shaft/err.no_connected_components.common.i

# A simulation of ActionComponents needs a mesh, which the classic simulation creates itself
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[ActionComponents]
  !include ../../components/shaft/err.no_connected_components.components.i
[]
