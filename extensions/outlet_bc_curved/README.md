# Curved Outlet Boundary Conditions

Copied and generalized SPARTA `surf_collide` styles for outlet boundaries on
arbitrary curved surfaces.

Included styles:

- `outletv/curved`
- `outletvparabolic/curved`
- `outletvgrad/curved`
- `outletvgradyfree/curved`

## Parameters

### `outletv/curved`

```text
vx vy vz
```

### `outletvgrad/curved`

```text
directionv directiongrad grad refvel width maxwidth
```

- `directionv`: velocity direction, `x=0`, `y=1`, `z=2`
- `directiongrad`: direction where the gradient is applied
- `grad`: velocity gradient
- `refvel`: reference velocity at the zero point
- `width`: width of the region with a gradient
- `maxwidth`: total system width

### `outletvparabolic/curved`

```text
meanv directionv maxx
```
