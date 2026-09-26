# Staging notes — read before pushing

This file lists judgment calls and open caveats from assembling this
repository. Nothing here has been pushed anywhere; review and resolve what
matters before `git push`.

## 1. `immunoscope/` commit provenance is unverified

`analysis.json` for each complex records `ImmunoScope 0.1.0, commit cbb0ae2b`
as the provenance of the analysis toolkit. The only copy of this code found
on disk is at `/data/work/development/immunoscope_recovered/immunoscope/`
(hence its directory name — it was apparently recovered, not a live clone),
and **it has no `.git` directory** — there is no way to confirm from this
copy alone that it is byte-identical to commit `cbb0ae2b`, or even that it
is the same version referenced by the analysis.json records used for this
manuscript's data.

Before citing this repository's `immunoscope/` as "the exact code used for
ITrace version 1.0," you should either:
- locate the original git history for ImmunoScope (if it exists anywhere
  else, e.g. another machine, an old remote, a different recovered path) and
  confirm/restore commit `cbb0ae2b`, or
- if this recovered copy is confirmed correct by other means (e.g. you
  recall recovering it directly from that state), say so explicitly in the
  manuscript rather than asserting a commit hash that can't be verified from
  this repo's own history, or
- if it's not the same version, decide what the correct citable version is.

Other candidate locations existed on disk but were not code copies (checked,
ruled out): `/srv/smbshare/manuscript/ImmunoScope` and `/srv/smbshare/doc/ImmunoScope`
are manuscript-drafting/report documents, not source. `../Immunex`, the
sibling-repo path named in `environment-imscope.yml`'s install instructions,
does not exist on this machine.

## 2. `pipeline/` is staged from an uncommitted working tree

`/data/work/Immuno-Dyn`'s git history stops at 2026-06-16; the copy of
`pipeline/` in this staging area is instead the **current, uncommitted
working-tree state** as of 2026-09-26, including 16 files that were never
tracked in that repo's git history at all (`essential_dynamics.py`,
`chain_roles.py`, `peptide_dihedrals.py`, `interface_interactions.py`,
`consensus.py`, `coupled_states.py`, `antigen.py`, `backfill_peptide_hla.py`,
`backfill_trust.py`, `fix_rmsf_regions.py`, `fix_run3_chain_mapping.py`,
`fix_run3_identity_metadata.py`, `tcr_cdr3.py`, `trust.py`, plus one renamed
file `recognition_mode.py` → `interface_descriptors.py`).

This staged copy is the actual code as it exists now, which is what you
asked for — but it means Immuno-Dyn's own git history no longer reflects
reality. Independent of this repo, you may want to commit that working-tree
state in `Immuno-Dyn` itself so your own project has a matching record.

Excluded from `pipeline/`: `run_analysis.py.bak_20260917_161438` (the only
backup file present; unambiguous).

## 3. `figs/` file selection

Included only the top-level analysis modules, `make_table*.py`, and the live
(non-`.bak`, non-archived-variant) scripts under `fig3/` and `fig4/`, per
the task scope (Technical Validation reproduction for Figures 3–4 and
Table 1 + supplementary tables). Excluded `fig1/`, `fig2/`, and everything
under `fig4/panels/`, `fig4/chimerax/`, `fig4/_v2_archive/`,
`fig4/_v3_candidate/`, `fig4/_v4_candidate/`, `fig3/_bak_20260903/`, all
`__pycache__/`, and all generated data/figure artifacts (`.tsv`, `.json`,
`.csv`, `.png`, `.pdf`, `.svg`) — those are outputs, not code, and are
already published as part of the manuscript's figures / the ITrace dataset
itself.

Spot-checked for stale/reverted files by comparing timestamps against
`.bak_*` siblings (`build_replicas.py`, `repdata.py`): in every case checked,
the live file is newer than all its backups, i.e. not a case of a `.bak`
representing a more current version than the live file. No ambiguous cases
found requiring your judgment.

If you also want Figure 1/2 (composition) or Table 1 lineage/provenance
figure-generation scripts included, they weren't copied — say so and they
can be added; `make_table1.py`, `make_table_lineage.py`,
`make_table_provenance.py`, and `make_table_reuse.py` themselves ARE already
included (top-level `figs/` scripts), only the Figure 1/2 plotting scripts
were excluded as out of scope for "Technical Validation reproduction."

## 4. Nothing else was excluded on uncertain judgment

Every file present in `pipeline/`, `immunoscope/`, and the selected `figs/`
scope was included as-is (minus `__pycache__`/`.pyc`/`.bak*`/generated
outputs, per the rules above). No files were dropped due to hesitation.
