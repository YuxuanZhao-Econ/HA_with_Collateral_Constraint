# Fixed forecast breaks and local equilibrium cache

`forecast_breaks_gh5.csv` contains the state-specific 18 breakpoints for the
20-coefficient cubic perceived law of motion. It is an input, not a generated
equilibrium. Keep it with the notebook/source when sharing the project.

The endpoints are the conditional minimum and maximum B, and the 16 interior
breaks are conditional quantiles r/17 (r=1,...,16), in the reference 12,000-date
simulation with burn=500 and seed=20260925. The source is the reference path
used by the completed `test/ha_spline_ks` experiment. `reference_n` counts that
reference path, not the final spline simulation. The knot positions stay fixed
throughout KS iteration. Columns `yT` and `yN` verify the state ordering.

`equilibrium_gh5_spline20_T48000.jls` is the notebook's reusable Julia Serialization
cache. It stores the model, fixed basis, coefficients, household solution,
simulation, iteration history, and a numerical-input signature. Generated
`.jls` files are ignored by Git. The old `equilibrium_gh5.jls` remains a record
of the previous nodal-rule experiment and is not loaded by the notebook.

The current notebook is configured for 48,000 training dates and an outer
stopping tolerance of 1e-3. It uses a separate filename from the earlier
12,000-date cache, so the earlier result is not silently reused. The fixed
breakpoints remain those of the original reference simulation.

The earlier `equilibrium_gh5_spline20.jls` cache was migrated from the completed 130-round, B-only,
12,000-date experiment after verifying an equivalent native Julia refit and
its stopping condition. It is not a new full KS execution. The notebook's
`FORCE_COMPUTE=false` loads the current 48,000-date cache when compatible;
`true` solves anew from the
affine initial rule. Changed model inputs, solver options, fixed breakpoints,
Julia version, dependencies, or numerical source invalidate the cache.
