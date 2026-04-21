# Curved Evaporation Boundary Conditions

Copied and generalized SPARTA `surf_collide` styles for evaporation and
condensation on arbitrary curved surfaces.

Included styles:

- `evapref/curved`
- `evaprefpart/curved`

## `evapref/curved`

```text
tsurf acc
```

- `tsurf`: surface temperature
- `acc`: condensation coefficient

Example:

```text
surf_collide inlet evapref 353.15 0.25
bound_modify ylo collide inlet
```

## `evaprefpart/curved`

Whole wall is treated as liquid:

```text
tsurf acc
```

Part of the wall is treated as liquid:

```text
tsurf acc direction liqmin liqmax
```

Example:

```text
variable coeff equal 0.25
variable nw equal ${coeff}*9.725e24

species ../air.species H2O
mixture water H2O vstream 0.0 0 0 temp 353.15 nrho ${nw}
region liquid block 0 5e-6 -1.5e-8 1.5e-8 INF INF
fix in emit/face water ylo region liquid
surf_collide inlet evaprefpart 353.15 ${coeff} 0 0.0 5e-6
bound_modify ylo collide inlet
```
