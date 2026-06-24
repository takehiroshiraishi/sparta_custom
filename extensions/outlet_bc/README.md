# Outlet Boundary Conditions

Custom SPARTA `surf_collide` styles for outlet boundaries.

Included styles:

- `outletv`
- `outletvparabolic`
- `outletvgrad`
- `outletvgradyfree`
- `outletvgradT`

## Parameters

### `outletv`

```text
vx vy vz
```

### `outletvgrad`

```text
directionv directiongrad grad refvel width maxwidth
```

- `directionv`: velocity direction, `x=0`, `y=1`, `z=2`
- `directiongrad`: direction where the gradient is applied
- `grad`: velocity gradient
- `refvel`: reference velocity at the zero point
- `width`: width of the region with a gradient
- `maxwidth`: total system width

### `outletvparabolic`

```text
meanv directionv maxx
```

### `outletvgradT`

Branching temperature-gradient outlet for planar box-face outlets, intended
first for 1D evaporation at `xhi`.

```text
t1 T1 t2 T2 uinf UINF [ut UT WT] [min_un V] [max_branch N] [tangent scale|resample] [rot scale|sample]
t1 T1 delta DELTAT uinf UINF [same optional keywords]
t1 T1 grad GRADT dx DXGHOST uinf UINF [same optional keywords]
```

- `T1`: scalar local interior outlet temperature for the current run.
- `T2`: scalar outside/ghost temperature. `delta` uses `T2 = T1 - DELTAT`; `grad`
  uses `T2 = T1 + GRADT*DXGHOST`.
- `UINF`: fixed outward bulk normal velocity.
- The first validated range is cooling only: `0 < T2 <= T1`.
- `max_branch` defaults to `20`; `min_un` defaults to `1.0e-12`.

Diagnostics are available through `stats_style sc_ID[N]`:

1. outgoing particles processed
2. incoming particles generated
3. sum of expected branch factors
4. mean expected branch factor
5. maximum expected branch factor
6. `uprime >= 0` count
7. `S2 < 0` count
8. invalid temperature count
9. branch cap count
10. `min_un` deletion count
11. current `T1`
12. current `T2`
13. safe-`U` diagnostic threshold
