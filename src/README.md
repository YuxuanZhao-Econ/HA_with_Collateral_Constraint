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

The HA implementation separates economic and numerical work from presentation:

- `model.jl` and `shocks25.jl`: primitives and grids.
- `forecasts.jl`: 20-coefficient cubic H, tangent tails, implied prices, and constrained regression.
- `households.jl`: household decisions and marginal values.
- `equilibrium.jl`: market clearing, distribution simulation, and the outer iteration.
- `impulse_responses.jl`: paired GIRF paths, conditional means, and CSV export.
- `plots.jl`: figures and table display, using explicit arguments.
- `cache.jl`: compatible equilibrium loading and saving.
- `validation.jl`: optional diagnostics, separate from the notebook's presentation.

Create additional files only when they contain an implementation. Add a
`common/` directory only if multiple models need shared numerical tools.

The entry file should load its supporting files relative to `@__DIR__`.
The notebook should load that entry file relative to the activated project root.

Local checks and exploratory variants belong in the Git-ignored `test/` or
`experiments/` directories. Reusable implementations needed to run a public
notebook belong in `src/`. Small exported summaries go in `results/`. The HA
notebook's reusable equilibrium cache goes in `Resources/ha_collateral_constraint/`,
where generated `.jls` files are ignored. Temporary checks remain in `tmp/`.
