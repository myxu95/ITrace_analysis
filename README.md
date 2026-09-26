# ITrace_analysis

Archived analysis-pipeline code for ITrace version 1.0, accompanying the
Scientific Data Data Descriptor submission for the ITrace pMHC–TCR molecular
dynamics trajectory resource. This repository covers the scientific analysis
pipeline only. The website and API implementation are maintained separately at
[github.com/myxu95/ITrace-web](https://github.com/myxu95/ITrace-web)
(MIT License).

## Structure

- **`protocol/`** — the GROMACS structure-preparation and equilibration
  protocol (`reproduce.sh` + `em.mdp`/`nvt.mdp`/`npt.mdp`/`md.mdp`): topology
  generation, solvation/ionization, energy minimization, NVT/NPT equilibration,
  production MD, and the trajectory-standardization post-processing
  (`gmx trjconv -pbc whole` → `-pbc nojump` → `-fit rot+trans`, protein-only
  extraction, topology-PDB dump).

- **`md_analysis/`** — the core per-system MD analysis toolkit: descriptor
  calculation (interface metrics, docking/incident angles, dihedral PCA, RMSF,
  essential dynamics, etc.). Version 0.1.0. This is the toolkit recorded in
  each complex's `analysis.json` under `provenance` as `md_analysis 0.1.0`.

- **`pipeline/`** — extraction, quality-control reporting, manifest building,
  and aggregation of `md_analysis`'s per-system outputs into the web-facing
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
- `environment-md-analysis.yml` (env `md-analysis`) — for the heavy per-system
  MD analysis (`md_analysis/`, driven via `pipeline/run_analysis.py`) and the
  publication figures (`figs/`).

`md_analysis` itself is not published to PyPI/conda; install it editable
inside the `md-analysis` environment:

```bash
conda env create -f environment-md-analysis.yml
conda activate md-analysis
pip install -e ./md_analysis
```

## Provenance

The per-complex `analysis.json` records provenance as `md_analysis 0.1.0`.
`md_analysis/` in this repository is the analysis code used to compute the
published per-system descriptors; no independently verifiable commit hash is
recorded for the exact state of that code, since its original version-control
history was not preserved separately from the working copy archived here.
`pipeline/` and `figs/` are tracked in the authors' own `Immuno-Dyn` git
history, unaffected by this caveat.

## License

MIT — see `LICENSE`.
