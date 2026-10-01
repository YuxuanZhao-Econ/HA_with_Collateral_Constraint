# Heterogeneous-agent collateral-constraint model

The model, calibration, numerical algorithm, and figures are presented in
`notebooks/HA_Collateral_Constraint.ipynb`. The notebook uses a 5×5 joint
aggregate-endowment chain, a shared 300-node individual/distribution asset grid,
and 160 mean-asset nodes. The notebook contains the explanations, explicit
experiment settings, calls, and displayed results; function definitions live here.

| File | Purpose |
| --- | --- |
| `shocks25.jl` | Five-point Gaussian-Hermite discretization in each aggregate dimension |
| `model.jl` | Parameters, grids, IID income, and consumption |
| `forecasts.jl` | Fixed cubic forecast, tangent tails, and constrained coefficient regression |
| `households.jl` | Household Euler time iteration and policies |
| `equilibrium.jl` | Nontradable market clearing, distribution simulation, and KS iteration |
| `steady_state.jl` | Constant-endowment household solution, invariant distribution, market-clearing price, and local cache |
| `validation.jl` | Household, market, resource, and independent Euler diagnostics |
| `impulse_responses.jl` | Paired generalized impulse responses, Monte Carlo errors, and CSV export |
| `plots.jl` | Steady-state, simulation, household-policy, aggregate-rule, and GIRF figures; table display |
| `cache.jl` | Load a compatible equilibrium or solve and save it |

The aggregate forecast is next-period mean assets, `H(z,B)`. Each state has
20 cubic B-spline coefficients: 17 polynomial intervals with matching values,
first derivatives, and second derivatives at their interior joins. Outside
the fixed spline interval, the endpoint tangent lines continue the rule.
The forecast is evaluated directly; the 160 B nodes are for household tables.
At every evaluation, the nontradable-price forecast follows from the
tradable resource constraint and nontradable market clearing:

```text
P(z,B) = (1-omega)/omega *
         ((R*B + yT[z] - H(z,B))/yN[z])^(1+eta).
```

The solver retains the full asset distribution during simulation. It clears
the current nontradable market at each date, regresses only the next mean
assets on visited states, and updates the 20 coefficients of `H` per state.
All fitting runs in Julia, without Python subprocesses or CSV exchange.
The notebook states
the grid, stopping condition, and interpretation of the resulting forecasts.

The fixed state-specific breakpoints are read from
`Resources/ha_collateral_constraint/forecast_breaks_gh5.csv`. They are the
conditional min/max and 16 conditional quantiles r/17 of a reference
12,000-date simulation, held fixed across outer iterations. A coefficient fit
minimizes squared B-next errors plus 0.002 times the squared change of the
forecast at the 160 household B nodes. It enforces aggregate-domain and
household-feasibility bounds over all B, including interior cubic extrema.
There is no node-index second-difference penalty or post-fit clipping.

## Notebook workflow

1. Section 5 activates the environment and includes `HACollateralConstraint.jl`.
2. Section 6 constructs the shared grids and calls `load_or_solve_steady_state`.
   It fixes both aggregate endowments at one, retains IID individual risk,
   and re-solves households without aggregate uncertainty. It displays
   stationary averages, equilibrium residuals, saving and consumption policies,
   and the invariant asset CDF. The result is cached in `steady_state.jls`;
   `FORCE_COMPUTE_SS=true` recomputes it.
3. Section 7 describes the stochastic numerical algorithm.
4. Section 8 calls `load_or_solve`. The current
   settings are 48,000 observations, a 500-date burn-in, and an outer stopping
   threshold of 1e-3, checked with tightly solved household policies.
5. Section 8.1 calls `simulate` for a separate 3,000-date path and
   `plot_simulation`, displaying the first 100 dates in the time-series panels.
6. Section 8.2 calls `plot_household_policies(result, z; B=B_policy)` for
   z=13 and z=12 at the same explicitly chosen B=-0.94. It interpolates the
   solved saving table and recovers consumption from the budget at the
   implied forecast price. These policy slices require no simulated distribution.
7. Section 8.3 calls `plot_aggregate_rule` for H and the implied log price.
   Panel titles count the final training-path observations; shading marks
   their conditional 1st–99th percentiles. Dotted vertical lines mark the fixed
   spline support endpoints; blue curves include the tangent tails.
8. Section 9 calls `generalized_irf`, `plot_girf`, and `save_girf_csv`.

The deterministic equilibrium is solved on individual assets and income only.
Household policy changes must be below 2e-10, invariant-distribution changes
below 1e-12, and the nontradable market residual below 1e-9. All averages and
the constrained share integrate over the stationary asset distribution and
the current IID income probabilities. It is distinct from the stochastic
economy's long-run average; the subsequent KS and GIRF initializations are
unchanged. Adding this separate solver does not invalidate a compatible
stochastic equilibrium cache.

Plotting functions return figures. The notebook chooses filenames, saves the
figures, and displays them. No plotting function reads notebook-global
`m`, `result`, `s`, or output-directory variables.

## Saved equilibria

`load_or_solve(m, cache_path; force_compute=false, solver_options...)` returns
the same result object as `solve_equilibrium`. The notebook uses
`Resources/ha_collateral_constraint/equilibrium_gh5_spline20_T48000.jls`, a local Julia
Serialization file. `FORCE_COMPUTE=true` computes a new solution; false reuses
a matching converged one and computes if none is available.

The cache signature includes the model arrays and parameters, supplied solver
settings, Julia version, numerical solver source, and Project/Manifest files.
The supplied solver settings include every fixed breakpoint. The new cache
format and filename distinguish it from the earlier nodal-rule equilibrium.
Plots, GIRFs, and optional diagnostics are excluded from the numerical-source
fingerprint, so editing presentation code does not force a new equilibrium solve.
An unconverged result is returned for inspection but is not saved as a reusable
equilibrium.

The earlier 12,000-date cache, `equilibrium_gh5_spline20.jls`, was converted from the completed 130-round experiment
in `test/ha_spline_ks/`. The new representation matches its H and log price
within 4e-13 on a dense grid; a native constrained refit gives the checked
update 5.9719e-4. This conversion does not rerun the full outer iteration.
The current 48,000-date configuration has a separate cache filename and does
not reuse that 12,000-date result. `FORCE_COMPUTE=true` starts afresh from the affine initial rule described in
the notebook, rather than from the experimental warm start.

## Generalized impulse responses

The two branches start from the same average inherited-asset distribution in
the final solver simulation. This is an unconditional time average across
aggregate states. The baseline impact state is z=13 and the adverse state is
z=12; default calculations use 40 dates and 2,000 paired paths.

Maximal coupling preserves each branch's Markov probabilities. A
fixed-distribution control variate reduces Monte Carlo noise. The four
endowment/price/consumption panels show percent differences of conditional
means; next mean assets are a level difference. The constrained-share panel
shows the two population percentages on a 0–100 scale. The CSV exports these
same quantities and their standard errors. The bands quantify simulation
sampling error, not equilibrium-approximation error.

## Accuracy and interpretation

The outer update and its two R-squared values are measured on the simulated
states. Both the iteration log and `forecast_diagnostics` use variation
within each aggregate state as the R-squared denominator; they report mean
assets and log prices separately.

Mean assets are an approximation to the aggregate distributional state.
Convergence and high forecast R-squared alone do not verify the exact
full-distribution equilibrium. Optional independent checks remain in
`validation.jl`; the current notebook does not execute them. Old diagnostic
and iteration CSVs in `results/` are not automatically refreshed.

The cubic basis imposes interior continuity and removes oscillatory cubic
extrapolation through tangent tails. It does not supply observations in
unvisited regions. The shaded visited interval does not guarantee dense
observations at every interior B, and rare states can remain weakly identified.
The forecast state remains aggregate shocks and mean assets; asset variance
is not an additional forecast state.
