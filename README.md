## Custom extension layout

Keep upstream SPARTA sources in `src/` as close to stock as possible.
Put local extensions under `custom/extensions/` and let the build stage
them into `src/` as symlinks before compiling.

Recommended layout:

```text
custom/
  README.md
  stage_extensions.sh
  extensions/
    outlet_bc/
      README.md
      src/
        surf_collide_outlet_v.cpp
        surf_collide_outlet_v.h
    my_solver/
      README.md
      src/
        collide_my_solver.cpp
        collide_my_solver.h
```

### Why this layout

- SPARTA's native build expects a flat `src/` directory.
- Style registration is filename-driven and scans headers in `src/`.
- GNU Make compiles `src/*.cpp`.
- CMake also globs `src/*.cpp` and `src/*.h`.

Because of that, custom code in nested directories will not be compiled
unless it is staged into `src/`.

### Naming rules

Use SPARTA's native filename and header registration conventions:

- `collide_<name>.cpp/.h`
- `surf_collide_<name>.cpp/.h`
- `react_<name>.cpp/.h`
- `compute_<name>.cpp/.h`
- `fix_<name>.cpp/.h`
- `region_<name>.cpp/.h`
- `dump_<name>.cpp/.h`
- command classes use the existing command naming pattern in `src/`

Each header still needs the usual registration block, for example:

```cpp
#ifdef SURF_COLLIDE_CLASS
SurfCollideStyle(mywall,SurfCollideMyWall)
#else
// class definition
#endif
```

### Build behavior

- `stage_extensions.sh` creates symlinks from `custom/extensions/**/src/*`
  into `src/`.
- `src/Makefile` should run the staging step before generating style headers.
- `src/CMakeLists.txt` should run the staging step during configure, before
  source globbing.

### Practical guidance

- Put extension sources under `custom/extensions/<package>/src/`.
- Keep package-local docs, examples, and notes next to the package, not in `src/`.
- Avoid editing upstream files unless you are changing a core interface.
- If you add or remove custom source files, rerun CMake configure before building.
- Do not define two different extension files with the same basename. Staging will
  fail on collisions.

### Existing legacy custom code

This tree already contains legacy custom code under `src/custom/...` with manual
symlinks into `src/`. That still works, but for new work prefer the top-level
`custom/extensions/` layout so custom code is not mixed into the upstream source
tree.
