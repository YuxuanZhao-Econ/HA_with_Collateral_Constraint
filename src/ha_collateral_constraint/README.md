# HA collateral constraint: mean-asset KS experiment

The reader-facing specification, algorithm, and executed results are in
`notebooks/HA_Collateral_Constraint.ipynb`. This is a decentralized-equilibrium
experiment with aggregate state `(z, B)`. The full asset histogram is retained
during simulation. No low-asset-share statistic is included in household forecasts.

| File | Purpose |
| --- | --- |
| `model.jl` | Bianchi parameters, IID income, grids, static CES allocation, forecast functions |
| `households.jl` | Lagged-policy Euler time iteration and current choices at a trial price |
| `equilibrium.jl` | Price clearing, Young lottery, fixed-path simulation, forecast fitting and damping |
| `coverage.jl` | Additional market-clearing observations for rare aggregate-state combinations |
| `experiment.jl` | Reproducible baseline, refinement, and coverage experiments with checked caches |
| `validation.jl` | Household KKT, market/distribution checks, recursive forecasts, independent equilibrium Euler checks |

## Economic inputs

Structural parameters and the full 16-state aggregate endowment chain are reused
from the Bianchi replication. The chain is read from its existing `shocks.jl`;
the original files are not edited. The new income shock is `epsilon = +/- 0.05`,
with equal probabilities, independent across households, time, and aggregate
shocks. This is an illustrative new parameter, not an estimated calibration.

## Iteration

For fixed price and mean-bond forecasting rules, iterate household policies using
the old saving choice inside expected marginal utility. At fixed current price,
CES homotheticity gives `lambda = K(p) * (cT)^(-sigma)`, so the scalar Euler
inversion is exact. The constrained branch and damping match the logic of the
replicated Bianchi solver. This is fixed-grid time iteration, not EGM.

Every simulated period recomputes household choices as the common nontradable
price is varied. Future marginal values and the forecast next mean remain fixed
during the current price-clearing calculation. The actual next distribution is
then computed from these choices. The forecast/actual mean discrepancy is recorded.

Forecast rules use a piecewise-linear basis in **signed** mean assets. Regressions
are separate across the 16 aggregate states. A small second-difference penalty
and prior on the preceding rule identify grid nodes outside the observed sample.
Forecast nodes are bounded to the numerical mean grid and to price ranges where
the entire individual policy grid admits positive expenditure and contains the
collateral limit. `forecast_guard_diagnostics` must be checked on validation paths.
These safeguards must not be mistaken for additional economic constraints.

Every outer round reuses the same aggregate draws and initial distribution.
Household policies warm-start from the preceding round, forecast updates are
damped, and household precision adapts to the previous rule-update norm. A
candidate outer stopping point is rechecked with tightly converged households.
Forecast convergence is measured on fitting states; changes across the entire
grid are reported separately. The notebook explains each step, the three asset
grids, all calibration values and their original targets, and the numerical
settings without requiring another notebook as a reference.

The third experiment adds aggregate-grid points near the collateral-sensitive
debt region and enriches the regression data. It selects actual simulated
histograms nearest the B nodes on its first round plus 100 regularly spaced dates.
Their date indices stay fixed across later rounds while their histograms evolve;
this avoids discontinuous changes in regression weights. It clears each
of their markets under every aggregate shock (all transition probabilities in
this calibration are positive). These unit-weight rows deliberately change the
projection measure toward rare states; they are not used for economic moments.
Its stopping criterion also includes these extra fitting states. Household
decisions still use only `(z, B)` for aggregate forecasting.

To reproduce the computation outside Jupyter, activate this repository's Julia
project, include `HACollateralConstraint.jl`, import its exports, and include
`experiment.jl`. Call `run_experiment(root; recompute=true)`, then
`run_coverage_experiment(root, runs[2]; recompute=true)`. The notebook also performs
the held-out validation and writes CSV diagnostics and figures.

## Accuracy

`converged` refers to the **projected forecast fixed point**, not an exact full-
distribution rational-expectations equilibrium. Report one-step and recursively
iterated forecast errors, held-out-shock-path errors, and independent Euler errors.
The latter clears next-period markets using the actual next histogram for each
possible aggregate shock, rather than simply reusing the expectation interpolant
that the household solver was designed to satisfy.

The notebook separates this approximation error from household KKT residuals,
market clearing, resource accounting, probability mass, and grid-boundary checks.
