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
  - heterogeneous households with IID individual income risk and aggregate shocks
  - Euler time iteration with collateral constraints
  - market clearing and simulation of the full asset distribution
  - aggregate forecasting rules based on the current shock and mean assets
  - baseline, refinement, and rare-state coverage experiments
  - held-out forecast errors and independent equilibrium Euler checks

- `src/`: model implementations (see [source layout](src/README.md))
  - `bianchi2011/Bianchi2011.jl`: representative-household replication
  - `ha_collateral_constraint/HACollateralConstraint.jl`: heterogeneous-household model
  - supporting files separate model inputs, household solvers, equilibrium,
    simulation, validation, and plotting

- `figures/`
  - exported notebook figures, organized by model

- `results/`
  - small CSV tables with simulation statistics, welfare, and numerical diagnostics
  - run metadata recording settings and source provenance

- `reference/`
  - Bianchi's article, online appendix, and supplied MATLAB replication package
  - the project's LaTeX model draft and compiled PDF

- `Project.toml` and `Manifest.toml`
  - direct dependencies and pinned transitive dependency versions

Local checks, exploratory experiments, and temporary solver caches belong in
`test/`, `experiments/`, and `tmp/`. These directories, including uppercase and
mixed-case variants, are excluded from Git and GitHub commits. Their contents
remain on disk and are not required in a fresh clone. The HA notebook recreates
missing caches and reuses them only when source and configuration checks match.

CSV files in `results/` preserve numerical tables for comparison and reuse
without rerunning the model. They are outputs, not solver inputs. The Bianchi
CSVs are earlier exported snapshots: its current notebook recomputes the
statistics on screen but does not refresh those files or their metadata.
The HA notebook exports its diagnostics, iteration history, and metadata when
its final results cell runs.

## Method

The Bianchi notebook uses the supplied structural parameters and 16-state joint
endowment process. DE and planner policies are computed by Euler time iteration:
use the preceding policies to form expected marginal values, check the
collateral constraint, update consumption and bonds, and recompute equilibrium
prices. The planner internalizes the effect of consumption on collateral
prices and updates the constraint multiplier as well.

The HA notebook follows a household-solution and aggregate-forecast iteration:

1. Guess the nontradable-price and next-period mean-asset forecasting rules.
2. Solve household consumption and saving policies given those rules.
3. Simulate aggregate shocks, clearing the nontradable market at each date.
4. Propagate the full asset distribution using the household choices.
5. Refit and damp the forecasting rules, then repeat until the projected rules converge.
6. Validate the resulting allocation on a held-out shock path, including
   independently re-cleared future markets for the Euler checks.

The notebooks distinguish policy-iteration convergence from equilibrium
accuracy. Budget constraints, collateral feasibility, market clearing,
distribution mass, numerical boundaries, and interpolation or forecasting
errors are checked separately.

## Interpretation

Aggregate uncertainty is represented by a Markov process. The Bianchi notebook
compares DE and planner outcomes under common shock paths, including the
transition costs of adopting the planner policy.

In the HA experiment, households forecast aggregate outcomes using the current
aggregate shock and mean assets. The simulation retains the full asset
distribution, but households do not condition their forecasts on that entire
distribution. Convergence of the projected forecasting rules therefore does
not establish an exact full-distribution rational-expectations equilibrium.
The current final candidate does not pass the notebook's 1% maximum-error
screen; the remaining approximation errors are reported explicitly.

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
