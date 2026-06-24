# Codex implementation spec: branching outlet temperature-gradient boundary in SPARTA

## Task

Implement a DSMC outlet boundary option in SPARTA that imposes an approximate temperature gradient at an outlet by converting outgoing molecules into incoming molecules using a shifted-Maxwellian mapping with stochastic branching.

This should be implemented as a new/custom outlet behavior rather than a conventional fixed-rate source. The method should not require explicit tracking of number density `n1`; it only requires a local outlet temperature estimate `T1`, a prescribed temperature drop/gradient, and a fixed bulk velocity `U_inf`.

## Physical assumptions

- Boundary normal points outward from the computational domain.
- A molecule crossing the outlet has outward normal velocity `un > 0`.
- The gas near the outlet is in local thermal equilibrium.
- The translational VDF is a shifted Maxwellian.
- The ideal gas law holds.
- Pressure is zero-gradient at the outlet.
- Temperature changes from inside value `T1` to outside/ghost value `T2`.
- Number density changes as `n2 = n1*T1/T2`, so `n1*T1 = n2*T2`.
- Bulk normal velocity is fixed as `U_inf` and should be in the outward-normal direction.

## Required user/input parameters

Use names consistent with the existing SPARTA code style, but the implementation needs these concepts:

- outlet boundary face/side, e.g. `xhi`, `yhi`, etc.
- prescribed temperature difference or gradient:
  - either `deltaT`, with `T2 = T1 - deltaT`, or
  - `gradT` and ghost distance `dx_ghost`, with `T2 = T1 + gradT*dx_ghost`.
- fixed outward bulk normal velocity `U_inf`.
- optional tangential bulk velocities `Ut`, `Wt`; default zero.
- local temperature source for `T1`: boundary cell translational temperature, preferably time-averaged.
- molecular specific gas constant `R = kB/m`, or equivalent in SPARTA internal units.
- optional numerical safeguards:
  - `min_T2 > 0`
  - `min_un` to avoid division by nearly zero
  - `max_branch` to prevent rare explosive branching
  - diagnostic counters.

Initial implementation may assume the common cooling case `0 < T2 < T1`. If `T2 >= T1`, either reject with a clear error/warning or add guarded behavior after validation.

## Core equations

For each outgoing particle, compute

```text
r = T2 / T1
A = 3 * R * T2 * log(T1 / T2)
S2 = r * (un - U_inf)^2 + A
uprime = U_inf - sqrt(S2)
```

where `un = dot(v, normal)` and `un > 0` for an outgoing particle.

The incoming branch is valid only if

```text
S2 >= 0
uprime < 0
```

Then define the expected branch multiplicity

```text
m = -uprime / un
```

Generate an integer number of incoming particles using stochastic rounding:

```text
N = floor(m)
if random_uniform_0_1() < (m - floor(m)):
    N += 1
```

This gives `E[N] = m`, so the incoming flux correction is correct in expectation.

## Remaining safe-U guideline

For `T2 < T1`, the colder outside Maxwellian can have a higher peak. To avoid the unmappable incoming central region, the input should satisfy approximately

```text
U_inf > sqrt(3 * R * T2 * log(T1 / T2))
```

For small `deltaT = T1 - T2`, this is approximately

```text
U_inf > sqrt(3 * R * deltaT)
```

The code should compute this diagnostic threshold and print/warn if `U_inf` is below it.

Do not enforce the older stricter condition

```text
U_inf > sqrt(3 * R * T2 * log(T1/T2) / (1 - T2/T1))
```

because branching removes the requirement `-uprime <= un`.

## Particle handling algorithm

When a particle crosses the selected outlet face:

1. Compute `un = dot(v, normal)`.
2. If `un <= 0`, do not apply this outlet conversion.
3. Get local/time-averaged outlet temperature `T1`.
4. Compute `T2` from the prescribed gradient or `deltaT`.
5. Guard:
   - if `T1 <= 0`, error;
   - if `T2 <= 0`, error;
   - if `un < min_un`, either delete the outgoing particle or use a guarded fallback;
   - if `S2 < 0`, delete outgoing particle and increment diagnostic counter.
6. Compute `uprime`.
7. Always remove/delete the original outgoing particle.
8. If `uprime >= 0`, generate no incoming particles.
9. If `uprime < 0`, compute `m = -uprime/un`.
10. Compute integer branch number `N` by stochastic rounding.
11. Optionally cap `N` by `max_branch`. If capped, increment a diagnostic counter because this introduces bias.
12. For `i = 1..N`, create one incoming particle just inside the boundary.

## Velocity assignment

Let the original velocity be decomposed as

```text
v = un*n + vt1*t1 + vt2*t2
```

where `n` is the outward normal and `t1,t2` are tangential basis vectors.

The new normal velocity is

```text
un_new = uprime     // negative means entering the domain
```

For deterministic tangential conversion, use

```text
vt1_new = Ut + sqrt(T2/T1) * (vt1 - Ut)
vt2_new = Wt + sqrt(T2/T1) * (vt2 - Wt)
```

If tangential bulk velocity is not used, set `Ut = Wt = 0`.

Reconstruct the full velocity:

```text
v_new = un_new*n + vt1_new*t1 + vt2_new*t2
```

Alternative optional mode: resample tangential components independently from a Gaussian at `T2`, using the existing SPARTA random thermal velocity utilities. This is more reservoir-like and reduces correlations, but deterministic scaling is simpler and closer to the current concept.

## Rotational energy for water

For water, use 3 rotational degrees of freedom. The simple deterministic mapping is

```text
erot_new = (T2/T1) * erot_old
```

This is unit-independent as long as `erot_old` and `erot_new` use the same SPARTA internal energy units.

For extra branched particles, deterministic copying/scaling produces correlated rotational energies. An optional better mode is to independently sample rotational energy at `T2` using the existing SPARTA rotational-energy sampling routine. For water with 3 rotational DOF, the equilibrium distribution is Gamma with shape `3/2` and scale `kB*T2` if energy is per molecule. Prefer existing SPARTA code paths for rotational energy units rather than hand-coding units.

## Position and timestep handling

Use the existing SPARTA pattern for boundary collisions or boundary particle insertion.

Recommended first implementation:

- place each incoming particle at the boundary crossing location, slightly inside the domain by a small epsilon along `-normal`;
- assign `v_new` immediately;
- follow the same remaining-timestep convention used by other SPARTA boundary/surface collision models.

If SPARTA already stores a crossing time or remaining fraction of the step, reuse that logic. Extra branched particles should use the same insertion position and remaining-time convention as the first converted particle.

## Diagnostics to add

Add counters and optional output for:

- number of outgoing particles processed;
- number of incoming particles generated;
- mean branch factor `m`;
- maximum branch factor observed;
- number of times `uprime >= 0`;
- number of times `S2 < 0`;
- number of times `T2 <= 0`;
- number of branch caps applied;
- current averaged `T1`, `T2`, and safe threshold `sqrt(3*R*T2*log(T1/T2))`.

These diagnostics are important for validating whether the branching boundary is stable and physically reasonable.

## Validation tests

Implement or run small tests before using in production:

1. `deltaT = 0`, `U_inf = 0`:
   - should reduce to specular-like normal reversal with unchanged tangential and rotational energy.

2. Small cooling gradient, e.g. `T2/T1 = 0.99`:
   - branch count should usually be near 0 or 1, with rare larger values.
   - outlet pressure should remain approximately flat.

3. Check sampled incoming distributions:
   - histogram of incoming normal velocity should be consistent with the intended mapped outside distribution over the accessible incoming range.
   - tangential temperatures should scale to `T2`.
   - rotational temperature should scale to `T2`.

4. Monitor macroscopic fields:
   - pressure near outlet should have weak/zero gradient;
   - temperature should develop the intended gradient;
   - density should vary approximately as `n ∝ 1/T` near the outlet.

## Important caveats

- This method avoids explicit `n1` tracking because the density ratio cancels analytically, but it still relies on actual outgoing molecules to drive incoming particle creation.
- The method is exact only in expectation, so it may increase statistical noise.
- Rare small-`un` outgoing particles can create large branch factors. Add diagnostics and optional caps.
- The method should first be implemented and validated for `T2 < T1`. Heating cases need separate checks because the square-root argument can become invalid for some outgoing velocities.
- The uploaded note contains a small typo: tangential scaling should be `sqrt(T2/T1)`, not `sqrt(T2/T2)`.
