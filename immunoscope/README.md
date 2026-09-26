# ImmunoScope Package Layout

This package is the product core of the repository.

## Layer Map

- `core/`: contracts, models, context, task discovery, exceptions
- `analysis/`: scientific and MD-domain logic
- `pipeline/`: workflow orchestration and batch execution
- `cli/`: user-facing command entrypoints
- `utils/`: compatibility helpers and generic support utilities only

## Boundary Rule

If a module defines product behavior, users should be able to reach it through this package without going through `scripts/` or `development/`.

## Dependency Direction

Preferred dependency direction is:

`cli/` -> `pipeline/` -> `analysis/` and `core/`


`analysis/` -> `core/` only when contracts or shared execution objects are needed

`utils/` should not grow into a second domain layer.

Historical implementations should stay outside the active product surface and
must not be reintroduced into `immunoscope/`.
