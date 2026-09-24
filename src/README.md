# Source Layout

Model directories use the lowercase notebook filename without `.ipynb`,
following the companion `sequence_space_jacobian` project.

| Notebook | Source directory | Entry file |
| --- | --- | --- |
| `Bianchi2011.ipynb` | `bianchi2011/` | `Bianchi2011.jl` |
| `HA_with_Collateral_Constraint.ipynb` | `ha_with_collateral_constraint/` | `HA_with_Collateral_Constraint.jl` |

The Bianchi baseline is implemented; see its [source audit and file map](bianchi2011/README.md).
The heterogeneous-household extension remains a reserved directory. Function definitions belong here;
the notebook contains the model explanation, algorithm explanation,
experiment settings, function calls, and displayed results.

As the implementation grows, split files by responsibility, for example:

- `parameters.jl`: parameter definitions.
- `grids.jl`: asset grids and income transitions.
- `household.jl`: household policies and borrowing constraints.
- `distribution.jl`: distribution updates and aggregation.
- `equilibrium.jl`: market clearing and equilibrium solution.
- `plots.jl`: reusable plotting functions.

Create these files only when they contain an implementation. Add a `common/`
directory only if multiple models need shared numerical tools. Numerical
methods and additional source files will be chosen when the model is implemented.

The entry file should load its supporting files relative to `@__DIR__`.
The notebook should load that entry file relative to the activated project root.
