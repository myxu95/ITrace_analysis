# Utils Package Policy

`immunoscope/utils` is a shrinking compatibility and helper layer.

Retired product interfaces should be removed from the active package tree
instead of being kept as a parallel product surface.

## Allowed Contents

- generic file/path helpers
- plotting helpers
- cleanup and selection helpers
- compatibility adapters that bridge old code to the new package layout

## Modules That Need Extra Care

Several current files in this package are domain-heavy and should eventually be reviewed for migration into clearer owners, especially chain identification and topology interpretation code.

## Ownership Clarifications

- CDR semantic detection is owned by the ANARCI-based topology analysis path, not by `utils/`

## Rule For New Code

Do not place new scientific or workflow-defining modules in `utils/` by default.

Before adding a file here, ask:

- is this generic enough to be shared by multiple subsystems?
- does it avoid owning scientific policy?
- would `analysis/`, `core/`, or `pipeline/` be a clearer owner?
