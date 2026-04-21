# Outlet Boundary Conditions

Custom SPARTA `surf_collide` styles for outlet boundaries.

Included styles:

- `outletv`
- `outletvparabolic`
- `outletvgrad`
- `outletvgradyfree`

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
