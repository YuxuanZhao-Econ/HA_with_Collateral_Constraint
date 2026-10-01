# Heterogeneous Agents with Collateral Constraints

This repository studies collateral-constrained borrowing in small open economies
through Julia notebooks. It begins with the representative-household model of
Bianchi (2011) and extends the framework to heterogeneous households with
individual income risk and aggregate endowment shocks.

The project is a teaching and research companion. Each notebook states the
economic environment and recursive equilibrium, explains its numerical method,
and presents policy functions, simulations, and economic and numerical checks.
Reusable numerical functions are kept in separate model directories.

## Requirements

The notebooks use Julia 1.12. The reproducible environment is recorded in
`Project.toml` and `Manifest.toml`; the recorded runs use Julia 1.12.4.
From the repository root, install the pinned dependencies with:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

Alternatively, from a Julia session:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

The dependencies include `IJulia`, `Plots`, `Roots`, and the standard libraries
listed in `Project.toml`. To open the project in Jupyter:

```julia
using IJulia
notebook(dir=pwd())
```

Open either notebook with the Julia 1.12 kernel and use **Run All**. Each notebook
activates the repository environment and loads its numerical helpers from
`src/` using project-root-relative paths. MATLAB is not required to run the
Julia notebooks.

## Project Structure

- [notebooks/Bianchi2011.ipynb](notebooks/Bianchi2011.ipynb)
  - representative-household small open economy with tradable and nontradable goods
  - decentralized equilibrium and constrained social planner
  - borrowing and consumption policies, collateral prices, and crisis dynamics
  - consumption-equivalent welfare gains and macroprudential instruments
  - equilibrium residuals, bond-grid refinement, and simulation sensitivity

- [notebooks/HA_Collateral_Constraint.ipynb](notebooks/HA_Collateral_Constraint.ipynb)
  - heterogeneous households with IID individual income risk and a 5×5 aggregate chain
  - deterministic aggregate steady state, stationary distribution, and household policies
  - Euler time iteration with collateral constraints
  - market clearing and simulation of the full asset distribution
  - a 20-coefficient cubic next-mean-asset rule H(z,B), with endpoint tangent tails and prices implied by aggregate resources
  - one configuration: 300 shared household/distribution nodes and 160 mean-asset nodes
  - simulated economy, household policies, aggregate rules, and generalized impulse responses

- `src/`: model implementations (see [source layout](src/README.md))
  - `bianchi2011/Bianchi2011.jl`: representative-household replication
  - `ha_collateral_constraint/HACollateralConstraint.jl`: heterogeneous-household model
  - supporting files separate model inputs, household solvers, equilibrium,
    simulation, validation, and plotting

- `figures/`
  - exported notebook figures, organized by model

- `results/`
  - small CSV tables with simulation statistics, welfare, and numerical diagnostics
  - the current HA GIRF export and historical experiment outputs

- `Resources/ha_collateral_constraint/`
  - fixed forecast breakpoints, `forecast_breaks_gh5.csv`, and separate deterministic (`steady_state.jls`) and stochastic (`equilibrium_gh5_spline20_T48000.jls`) Julia equilibrium caches
  - generated `.jls` files are kept locally and excluded from Git

- `reference/`
  - Bianchi's article, online appendix, and supplied MATLAB replication package
  - the project's LaTeX model draft and compiled PDF

- `Project.toml` and `Manifest.toml`
  - direct dependencies and pinned transitive dependency versions

Local checks, exploratory experiments, and temporary files belong in
`test/`, `experiments/`, and `tmp/`. These directories, including uppercase and
mixed-case variants, are excluded from Git and GitHub commits. Their contents
remain on disk and are not required in a fresh clone.

Section 6 of the HA notebook solves a deterministic aggregate steady state,
keeping both endowments at one while retaining IID individual income risk.
It reports stationary averages, an asset CDF, and household saving and consumption
policies. `FORCE_COMPUTE_SS=false` reuses its compatible local cache.
Section 7 introduces the stochastic algorithm.

In Section 8, `FORCE_COMPUTE=false` loads a compatible cached
equilibrium or computes and saves one if needed. Set it to `true` to recompute.
Compatibility checks cover the model, solver settings, Julia version, numerical
source, and pinned dependencies. Editing figures does not invalidate the
equilibrium cache. Only converged results are saved.

CSV files in `results/` preserve numerical tables for comparison and reuse
without rerunning the model. They are outputs, not solver inputs. The Bianchi
CSVs are earlier exported snapshots: its current notebook recomputes the
statistics on screen but does not refresh those files or their metadata.
The HA notebook currently exports the GIRF responses, the baseline and shocked
constrained-household shares, and their Monte Carlo standard errors. Existing
HA diagnostics and iteration-history CSVs are historical outputs; the current
notebook does not refresh them.

## Method

The Bianchi notebook uses the supplied structural parameters and 16-state joint
endowment process. DE and planner policies are computed by Euler time iteration:
use the preceding policies to form expected marginal values, check the
collateral constraint, update consumption and bonds, and recompute equilibrium
prices. The planner internalizes the effect of consumption on collateral
prices and updates the constraint multiplier as well.

After the deterministic reference equilibrium, the HA notebook follows a
household-solution and aggregate-forecast iteration:

1. Guess H(z,B), a fixed 20-coefficient cubic next-period mean-asset rule with
   tangent-line tails, and derive its price forecast
   from the tradable resource constraint and nontradable market clearing.
2. Solve household consumption and saving policies given H and its implied prices.
3. Simulate aggregate shocks, clearing the nontradable market at each date.
4. Propagate the full asset distribution using the household choices.
5. Refit and damp H using only observations on the simulated path. Repeat until
   the changes in H and its implied log price on those states are below 1e-3.
6. Display a separate simulation, household policy slices, fitted aggregate
   curves with state observation counts, and a generalized impulse response.

The HA fit is configured for 48,000 dates after a 500-date burn-in and prints separate
within-state R-squared values for next mean assets and log prices. Its GIRF
compares z=13 with z=12 from a common time-averaged asset distribution. The
constrained-share panel displays the two population percentages. The separate
3,000-date simulation shows its first 100 dates in the time-series panels.

Independent HA Euler and feasibility diagnostics remain available in
`src/ha_collateral_constraint/validation.jl`; the notebook does not run them.

## Interpretation

Aggregate uncertainty is represented by a Markov process. The Bianchi notebook
compares DE and planner outcomes under common shock paths, including the
transition costs of adopting the planner policy.

In the HA experiment, households forecast aggregate outcomes using the current
aggregate shock and mean assets. The simulation retains the full asset
distribution, but households do not condition their forecasts on that entire
distribution. Convergence of the projected forecasting rules therefore does
not establish an exact full-distribution rational-expectations equilibrium.
High forecast R-squared values and a converged iteration are not a substitute
for independent equilibrium and grid-accuracy checks.

## Scope and Extensions

- The Bianchi notebook is an equation-based quantitative replication. Its
  parameters and shock arrays come from the supplied MATLAB package, but
  numerical updates and simulation settings differ. The appendix describes
  value function iteration for the planner, while the Julia implementation
  iterates its first-order conditions. See the
  [source audit](src/bianchi2011/README.md) for the implementation differences.
- Actual MATLAB execution results have not been verified. Local Julia
  experiments isolating a change in the planner multiplier are not executions
  of the complete MATLAB script.
- The HA notebook currently solves a decentralized-equilibrium approximation.
  Individual income-risk parameters are illustrative. An HA planner and
  optimal macroprudential policy are outside the current implementation.
- Improving the HA forecast approximation and validating its accuracy remain
  necessary before drawing precise quantitative policy conclusions. See the
  [HA algorithm and accuracy notes](src/ha_collateral_constraint/README.md).

## References and Resources

- Bianchi, J. (2011). Overborrowing and Systemic Externalities in the Business
  Cycle. *American Economic Review*, 101(7), 3400-3426.
  [Article](https://doi.org/10.1257/aer.101.7.3400).
- Local article, online appendix, and original MATLAB replication materials:
  `reference/`.
