# ITrace_analysis

Archived analysis-pipeline code for ITrace version 1.0, accompanying the
Scientific Data Data Descriptor submission for the ITrace pMHC–TCR molecular
dynamics trajectory resource. This repository covers the scientific analysis
pipeline only. The website and API implementation are maintained separately at
[github.com/myxu95/ImmunoTrace-web](https://github.com/myxu95/ImmunoTrace-web)
(MIT License).

**Before pushing this repository, read `NOTES.md` — it lists provenance
caveats that need a decision.**

## Structure

- **`immunoscope/`** — the core per-system MD analysis toolkit: structure
  preparation, trajectory standardization, and descriptor calculation
  (interface metrics, docking/incident angles, dihedral PCA, RMSF, essential
  dynamics, etc.). Version 0.1.0. This is the toolkit recorded in each
  complex's `analysis.json` under `provenance` as `ImmunoScope 0.1.0`.

- **`pipeline/`** — extraction, quality-control reporting, manifest building,
  and aggregation of ImmunoScope's per-system outputs into the web-facing
  `analysis.json` records and `web_data/` bundle. Includes structural
  annotation (chain roles, TCR/antigen identity, CDR contacts), QC gating,
  and the post-hoc correction scripts applied during preparation of this
  release (e.g. `fix_run3_chain_mapping.py`, `fix_rmsf_regions.py`).
  `pipeline/tests/` holds the corresponding unit tests.

- **`figs/`** — the statistics and figure-generation scripts used to produce
  the manuscript's Technical Validation figures and tables (Figure 3, Figure
  4, Table 1, and the composition/reuse/lineage/provenance tables), reading
  from the aggregated per-complex data produced by `pipeline/`.

## Environment

Two conda environments are used:

- `environment.yml` (env `immuno-web`) — for `pipeline/` (extract / enrich /
  metrics / manifest / thumbnails / aggregate).
- `environment-imscope.yml` (env `imscope`) — for the heavy per-system MD
  analysis (`immunoscope/`, driven via `pipeline/run_analysis.py`) and the
  publication figures (`figs/`).

`immunoscope` itself is not published to PyPI/conda; install it editable
inside the `imscope` environment:

```bash
conda env create -f environment-imscope.yml
conda activate imscope
pip install -e ./immunoscope
```

## Provenance

The per-complex `analysis.json` records provenance as `ImmunoScope 0.1.0,
commit cbb0ae2b`. See `NOTES.md` for an open caveat on verifying that this
archived copy of `immunoscope/` matches that exact commit.

## License

MIT — see `LICENSE`.
