# Source Layout

Model directories use the lowercase notebook filename without `.ipynb`,
following the companion `sequence_space_jacobian` project.

| Notebook | Source directory | Entry file |
| --- | --- | --- |
| `Bianchi2011.ipynb` | `bianchi2011/` | `Bianchi2011.jl` |
| `HA_Collateral_Constraint.ipynb` | `ha_collateral_constraint/` | `HACollateralConstraint.jl` |

The Bianchi baseline is implemented; see its [source audit and file map](bianchi2011/README.md).
The heterogeneous-household DE experiment is also implemented; see its
[algorithm and accuracy notes](ha_collateral_constraint/README.md). It uses
aggregate-shock and mean-asset forecasts while propagating the full asset
distribution during simulation. Convergence of the projected forecasts does
not establish accuracy of the full-distribution equilibrium.

Function definitions belong here; the notebooks contain the model and algorithm
explanations, experiment settings, function calls, and displayed results.

Split files by responsibility, for example:

- `parameters.jl`: parameter definitions.
- `grids.jl`: asset grids and income transitions.
- `household.jl`: household policies and borrowing constraints.
- `distribution.jl`: distribution updates and aggregation.
- `equilibrium.jl`: market clearing and equilibrium solution.
- `plots.jl`: reusable plotting functions.

Create additional files only when they contain an implementation. Add a
`common/` directory only if multiple models need shared numerical tools.

The entry file should load its supporting files relative to `@__DIR__`.
The notebook should load that entry file relative to the activated project root.

Local checks and exploratory variants belong in the Git-ignored `test/` or
`experiments/` directories. Reusable implementations needed to run a public
notebook belong in `src/`. Small exported summaries go in `results/`; temporary
solver caches go in the ignored `tmp/` directory.
