# GapHeatTransfer

## Description

This BC applies heat transfer across unmeshed gaps between different blocks.

The `quadrature` option is generally recommended for most models. With this option, heat flux is
computed and applied as an integrated boundary condition at the quadrature points on both faces. Use
of the quadrature options generally gives smoother results, although there can be small differences
in the heat flux on the two surfaces.

It is also important to use the appropriate `gap_geometry_type` parameter (PLATE, CYLINDER, or
SPHERE) for the model geometry.

Two-dimensional Cartesian geometries are not restricted to be in or parallel to the X-Y coordinate plane.

The heat flux $q$ applied by this BC is computed as follows:

!equation
q = f h (T_B - T_A) \,,

where

- $T_A$ is the surface temperature of the current boundary,
- $T_B$ is the surface temperature of the coupled boundary,
- $h$ is the thermal conductance, and
- $f$ is the "edge multiplier".

This BC has two formulations:

- The "quadrature-based" formulation
- The "node-based" formulation

The formulation is selected via the [!param](/BCs/GapHeatTransfer/quadrature)
parameter (setting to `true` gives the "quadrature-based" formulation). The
formulation determines how to compute $T_B$ and $f$.

The thermal conductance $h$ is provided by a regular material property with the
name `gap_conductance<suffix>`, where `suffix` is the value of the parameter
[!param](/BCs/GapHeatTransfer/appended_property_name).

### Quadrature-based Formulation

For the quadrature-based formulation, the [PenetrationLocator.md] is used to
get information on the coupled boundary for the current quadrature point $\mathbf{x}_q$
(in real space). The coupled boundary temperature $T_B$ is then computed as

!equation
T_B(\mathbf{x}_q) = \sum\limits_{j\in\mathcal{I}} T_j \phi_j(\mathbf{x}_B(\mathbf{x}_q)) \,,

where $\mathbf{x}_B(\mathbf{x}_q)$ is the projected position of $\mathbf{x}_q$
on the coupled boundary and $\mathcal{I}$ is the set of degree of freedom indices
on the face containing $\mathbf{x}_B(\mathbf{x}_q)$.

### Node-based Formulation

For the node-based formulation, the gap temperature and distance are provided
via the parameters [!param](/BCs/GapHeatTransfer/gap_temp) and [!param](/BCs/GapHeatTransfer/gap_distance),
respectively. The edge multiplier is set to unity: $f = 1$.

## Example Input syntax

!listing modules/heat_transfer/test/tests/heat_conduction/2d_quadrature_gap_heat_transfer/nonmatching.i block=ThermalContact

!syntax parameters /BCs/GapHeatTransfer

!syntax inputs /BCs/GapHeatTransfer

!syntax children /BCs/GapHeatTransfer
