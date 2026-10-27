!include ../../components/component/err.nonexisting_component.common.i

# A simulation of ActionComponents needs a mesh, which the classic simulation creates itself
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[ActionComponents]
  !include ../../components/component/err.nonexisting_component.components.i
[]
