# Bianchi (2011) replication

`notebooks/Bianchi2011.ipynb` is the reader-facing model, algorithm, and results.
Load numerical functions through `Bianchi2011.jl`. Plot helpers are loaded
separately from `plots.jl` so numerical checks do not need a plotting backend.

| File | Responsibility |
| --- | --- |
| `model.jl` | Calibration, preferences, prices, bond grid, interpolation |
| `shocks.jl` | Exact Float64 export of the supplied 16-state Markov process |
| `export_shocks.py` | One-time SciPy conversion; not required by the notebook |
| `solvers.jl` | DE and planner Euler time iteration; fixed-policy welfare evaluation |
| `simulation.jl` | Common-shock simulation, crisis definition, matched crisis counterfactual, welfare and tax |
| `validation.jl` | Grid and off-grid residuals, derivative checks, cross-grid policy distances |
| `plots.jl` | Reusable figure generation and PNG/SVG export |

## Source audit and intentional differences

Sources are the supplied AER article, online appendix, and MATLAB package under
`reference/`. Original reference files are not edited. The notebook verifies
the SHA-256 of `proc_shock.mat` before using its Julia export.

- **Calibration:** use `main.m:23-28` (beta 0.906, omega 0.307, kappa 0.3235),
  not the paper's rounded table values. The joint shock process is reused,
  not re-estimated from data or replaced by independent shocks.
- **Grid:** the MATLAB script sets 100 points (`main.m:63`); Appendix section 2
  specifies 800. The main notebook uses 800, with 400/1600-point checks on the
  same [-1.02, -0.2] interval. Neither bond-grid endpoint is reached in the
  baseline simulations.
- **DE:** Euler time iteration with expected marginal utility evaluated at
  the previous bond policy, as in `main.m:138`. Each iteration calculates
  the borrowing limit algebraically from the previous price. The binding
  check precedes scalar Euler inversion with `Roots.find_zero(..., Brent())`.
  Consumption/bonds are damped with
  weight 0.2, then prices are recomputed. There is no precomputed collateral
  boundary and no candidate-dependent expectation inside the scalar root.
- **Planner:** the same iteration with marginal value `uT + mu*Psi` and
  consumption/bond damping 0.005 (the script uses 0.02 on its 100-point grid). At a binding choice, update the multiplier
  as `max(0, lambda_old - expected)`; at an interior choice it is zero.
  The multiplier update is retained without policy damping, as in the script.
  This implements the model's Euler and complementary-slackness conditions.
  The extra beta in `main.m:260` (the expectation at line 256 already includes
  beta) and the stale root-start variable at line 277 are not reproduced.
- **Update convention:** after damping consumption, recover bonds from the
  budget and prices from market clearing. This keeps those identities exact
  at each iteration; unlike the script, there is no additional undamped
  collateral projection. The collateral constraint is enforced at convergence.
  Stop on undamped policy/multiplier updates below 1e-11 and independently
  check the final equilibrium residuals. The MATLAB stopping tolerance is 1e-5.
- **Interpolation:** use the converged bond policy's linear interpolant in
  simulation, with a 1e-5 collateral-slack tolerance for classifying binding
  observations. No off-grid optimization or boundary projection is performed.
  Report off-grid Euler errors and collateral slack, and refine the grid.
- **Welfare:** evaluate the final policies with fresh flow utility arrays.
  The MATLAB variables `totalcDE` and `totalcSP` were last computed inside
  their solver loops, before the final update. Value evaluation includes
  the transition from a common initial state to the planner allocation.
- **Crises:** use a binding constraint AND net outflows above DE mean plus
  one DE standard deviation. This agrees with the paper and executed code
  (`main.m:519`), rather than that line's misleading "2 sd" comment.
  Date t binding is matched to date t's choice b[t+1]; no shifted date test.
- **Table 2:** select the median DE crisis by raw net outflow b[t+1]-b[t].
  Restart the planner two periods before that event at the same DE debt and
  apply the same shocks, as in footnote 12. The visible MATLAB summary instead
  averages SP observations from its separate history at DE crisis dates.
  Single-event numbers depend on the selected simulated episode.
- **Simulation:** use a documented Xoshiro seed, 5,000 burn-in periods and
  80,000 retained observations. The supplied script uses 79,000 observations
  with no burn-in and does not set an RNG seed. Its commented reference to
  `shocks_initials_aer.mat` points to a file absent from the supplied package.
- **Real exchange rate:** report the domestic CES price index in tradable
  units, as computed in `main.m:496-497`; its decline is reported as positive
  depreciation. The inverse exponent printed in the article's footnote 13
  conflicts with its accompanying positive-comovement statement and code.
  This convention is explicit in the notebook and all figures.
- **Densities:** histograms integrate to one, conditional on the relevant
  sample. The original plotting script divides KDE values by simulation
  length; that vertical normalization is not reproduced.

This is an equation-based quantitative replication, not a claim of bitwise
agreement with an execution of the unmodified MATLAB package. Remaining
differences from published targets are shown rather than recalibrated away.
