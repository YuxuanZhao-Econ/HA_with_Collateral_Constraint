# Heterogeneous Agents with Collateral Constraints

This project studies a heterogeneous-household extension of Bianchi (2011),
*Overborrowing and Systemic Externalities in the Business Cycle*. The intended
environment is a small open economy with tradable and nontradable goods,
household income risk, and borrowing limits tied to the market value of income.

The first research notebook, [Bianchi2011.ipynb](notebooks/Bianchi2011.ipynb),
replicates the original representative-household model before introducing
heterogeneity. It defines the model and equilibrium, explains the algorithms,
and compares policies, crises, welfare, and macroprudential instruments with
the published results.

## Requirements

The project uses Julia 1.12. `Project.toml` and `Manifest.toml` record the
environment based on the companion `sequence_space_jacobian` project, with
the replication's standard-library dependencies added. The manifest uses
Julia 1.12.4. Direct dependencies are `IJulia`, `Plots`, `LinearAlgebra`,
`Printf`, `Random`, `Statistics`, and `SHA`.

From the repository root, install the pinned dependencies with:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

Alternatively, from a Julia session started in the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using IJulia
notebook(dir=pwd())
```

## Project Structure

```text
HA_with_Collateral_Constraint/
|-- README.md
|-- Project.toml
|-- Manifest.toml
|-- .gitignore
|-- notebooks/
|   `-- Bianchi2011.ipynb         Original-model replication and exposition
|-- src/
|   |-- README.md               Source organization conventions
|   |-- bianchi2011/             Model, solvers, simulation, checks, and plots
|   `-- ha_with_collateral_constraint/   Reserved for the HA extension
|-- figures/bianchi2011/         Notebook-generated PNG and SVG figures
|-- results/bianchi2011/         Comparison tables and execution provenance
|-- reference/                  Papers and original replication materials
|-- main.tex                    Existing model draft
|-- Trade_Final_Project.pdf     Existing compiled draft
|-- experiments/               Local experiments; ignored by Git
|   `-- consumption_saving/     EGM, borrowing limits, and asset-grid boundary experiments
|-- test/                       Local numerical checks; ignored by Git
`-- tmp/                        Temporary files and logs; ignored by Git
```

The subsequent HA notebook, `notebooks/HA_with_Collateral_Constraint.ipynb`,
remains reserved. Its functions will live in the corresponding lowercase
source directory. Local ignored directories can be recreated after cloning:

```julia
mkpath.(joinpath.(pwd(), ["experiments", "test", "tmp"]))
```

The local [consumption-saving experiments](experiments/consumption_saving/README.md)
collect the fixed-interest-rate household exercises, their
[discussion summary](experiments/consumption_saving/DISCUSSION.md), Julia scripts,
figures, and numerical checks. These are partial-equilibrium household problems,
not an RBC model or the Bianchi replication. This ignored directory is available
locally and is not included in a Git clone.

## Notebook and Source Workflow

The main notebook first describes the model, including timing, household
states and choices, income processes, collateral constraints, and equilibrium.
It then explains the numerical algorithm before presenting computations
and diagnostics. Calibration choices and research exercises remain visible in
the notebook.

Reusable function definitions belong in `src/`. Following the companion
project, the model directory uses the lowercase notebook filename without the
`.ipynb` extension. The replication entry file is
`src/bianchi2011/Bianchi2011.jl`. See [source layout](src/README.md).

The notebook locates and activates the repository's `Project.toml`, then
loads source files using paths relative to that root. Experiments, checks, and
temporary outputs belong in the corresponding local directories rather than
the main notebook's source directory.

## Replication Method and Scope

The notebook uses the unrounded calibration and 16-state joint endowment
process from the supplied MATLAB package. Its main bond grid has 800 points,
as in the online appendix. Both competitive equilibrium and the constrained
planner use Euler time iteration: each iteration fixes expected marginal
values at the old bond choice, calculates the borrowing limit from the old
price, chooses consumption, and updates policies and prices. The planner also
updates the collateral multiplier. Future aggregate shocks enter all
conditional expectations.

The notebook reports:

- Borrowing policies, collateral prices, and the ergodic debt distribution.
- Crisis frequencies and unconditional moments compared with the paper.
- The matched median-crisis counterfactual from footnote 12.
- Consumption-equivalent welfare gains and implied debt taxes.
- Grid and off-grid residuals, 400/800/1600-point comparisons, and simulations
  with multiple fixed seeds.

The default simulation keeps 80,000 observations after a 5,000-period burn-in.
Longer robustness simulations keep 500,000 observations for each seed.
Figures and small result tables are regenerated by running all notebook cells.
Timing depends on compilation and hardware; the 1,600-point planner check
allocates roughly 0.7 GB for its utility and marginal-utility tables.

This is an equation-based quantitative replication, not a bitwise MATLAB port.
The [source audit](src/bianchi2011/README.md) records discrepancies between the
paper, appendix, and distributed script, including the planner loop and crisis
event comparison. Remaining published-target gaps are reported without
recalibration. The fixed-tax exercise and parameter-sensitivity table are
outside this baseline. Household heterogeneity is the next stage.

Open the notebook with the Julia 1.12 kernel and use **Run All**. Optionally,
Python with `nbconvert` and `nbclient` can execute it from the repository root:

```bash
python -m jupyter nbconvert --execute --to notebook --inplace --ExecutePreprocessor.timeout=1800 notebooks/Bianchi2011.ipynb
```

Notebook execution does not require SciPy or MATLAB: the original shock file
has a full-precision Julia export and its source hash is checked at runtime.

## Version Control

Track the notebook with executed outputs, source files, documentation, selected
figures, small result tables, and both Julia environment files. The `.gitignore` excludes local experiment,
test, and temporary directories, including singular/plural and uppercase
variants, as well as notebook checkpoints and generated runtime files.

The existing LaTeX draft, compiled PDF, reference papers, and original MATLAB
replication package remain in their current locations.

## References

- Bianchi, J. (2011). Overborrowing and Systemic Externalities in the Business
  Cycle. *American Economic Review*, 101(7), 3400-3426.
  [Article](https://doi.org/10.1257/aer.101.7.3400).
- Local paper, online appendix, and MATLAB replication code: `reference/`.
