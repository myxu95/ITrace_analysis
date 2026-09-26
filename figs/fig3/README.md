# Figure 3 — what the trajectories say about the peptide

Two figures, each drawn on its own at 141 mm (Scientific Data full text width)
and each sized by its own content.  They replace panels b and c of the old NAR
`figures_nar/fig4_crosssystem.png`, redrawn in house style against the current
release set; how they are combined into the printed Figure 3 is a separate
decision.

Neither carries descriptive prose.  A journal figure is an image and its legend
is manuscript text, so every sentence lives in `../FIGURE_LEGENDS.md` and only
what is needed to READ the graph is drawn: axis titles, position labels, counts
and the one data strip under each plot.

| File | Size | Shows |
|---|---|---|
| `fig_peptide_position` | 141 × 99 mm | TCR contact occupancy and residue SASA at P1–P9, 154 nine-mers, plus the anchor-flag rate per position |
| `fig_bulge_by_length` | 141 × 72 mm | peptide bulge height against peptide length, all 245 complexes |

Counting unit is the **complex** (245), not the trajectory (735): the three
replicas of a complex start from the same structure, so a replica is not an
independent observation of the library.  Every value plotted is already the
mean over a complex's replicas, and every box therefore shows spread ACROSS
COMPLEXES.

## What changed from the NAR panels

- The old panel b put mean TCR contact on a left axis and buried SASA on a
  reversed right axis labelled "(rev.)", with both axis titles set in the
  colour of their series.  Both quantities are strongly bimodal across the
  library, so a mean landed where almost no complex sits; they are now two
  stacked box rows on one shared P1–P9 axis, with achromatic axis titles and no
  inverted axis anywhere.
- MHC contact occupancy is not plotted at all.  It saturates (median
  0.95–1.00 at every position), which is why the burial axis is SASA.
- The `anchor` flag is no longer implied by the shape of two curves.  It is
  printed as a number per position, and its two fixed cut-offs are named in the
  legend — they are absolute thresholds on measured quantities, not quantiles
  of this library.
- The old panel c dropped any peptide length represented by fewer than five
  complexes, which silently discarded the 4-mer and the 12-mer.  Both are now
  drawn as the single complex they are and named, so the figure accounts for
  all 245.
- The axis break between 4 and 8 residues matches the retired
  `../_retired/fig_peptide_length.py`, drawn the same way for the same reason.
- The source is `web_data_1000`, the post-dedup release set.  The NAR script
  still points at the superseded `web_data`.

## Rebuild

    python ../build_dynamics.py       # only when the release set changes
    python fig_peptide_position.py
    python fig_bulge_by_length.py
    python ../build_legends.py        # after either of the above

Each writes `<name>.pdf` (vector, for submission) and `<name>.png` (600 dpi
RGB, for reading) next to itself, then prints the rendered size and an audit of
any text that left the canvas.  That audit line must stay empty.

## Where the inputs and the style live

One level up, shared with every other figure folder:

- `../figstyle.py` — canvas, type scale, colours, `header()`, `box()`,
  `render()`, `save()`
- `../dyndata.py` — reads the two frozen tables below, and owns the single
  quantile definition every box in this folder uses
- `../dynamics_bulge.tsv` — 245 rows, one per complex: peptide length and
  replica-mean bulge height
- `../dynamics_pep9.tsv` — 1386 rows, one per (9-mer complex, position):
  replica-mean `tcr_contact`, `hla_contact` and `sasa_nm2`

Both tables are built by `../build_dynamics.py` out of
`/home/xmy/work/data/immunotrace/web_data_1000`, reading `analysis.json`
`geometry.bulge_height_angstrom` and `interface.peptide_table`.  Freezing them
means the figures rebuild without the 23 GB release mount and every number in
Figure 3 stays auditable from the repository alone.

## Legends

Each script keeps its prose in a module-level `LEGEND` dict — title, subtitle
and notes — built inside `draw()` because it is full of computed numbers.
`figstyle.LEAN` (default `True`) keeps it off the canvas; `../build_legends.py`
imports every script, harvests the same dict and writes `../FIGURE_LEGENDS.md`.
Edit the wording in the script, never in the markdown, and re-run
`build_legends.py`.

Set `figstyle.LEAN = False` to render an annotated copy for internal review.
The claims those notes make are guarded where they can be: `fig_bulge_by_length`
asserts that the 8-, 9-, 10- and 11-residue interquartile ranges really are
disjoint before the note says so.

## House rules these figures follow

Nothing is set below 6.5 pt and nothing is set in grey.  Colour encodes data:
`PRIMARY` for the library distribution, and `SECOND` for nothing but the
exception class -- 7byd and 6zkx, the lengths represented by a single complex,
in the same ochre 7byd was given in the retired `_retired/fig_peptide_length.py`,
and the four reverse-polarity complexes in `anglestrip.py`.  The lower row of
`fig_peptide_position` used to take that ochre too, which was wrong twice over:
solvent accessibility is a second measurement of the whole library and not an
exception to it, and on the composed sheet it landed a panel away from the
reverse-polarity ochre and read as the same class.  It is now `THIRD` -- the
teal added to `figstyle` for exactly this case, Fig 1's contour hue at the
weight `PRIMARY` and `SECOND` share -- so the sheet carries one hue per role:
blue for the population, teal for its second measurement, ochre for anything
exceptional.  Full spec in `../../图风格统一规范.md`.

## Combining them later

Both scripts expose `draw(ax, top, w)`: they lay out downwards from `top` and
return the y of their lowest ink.  A composed sheet is a new script that calls
them on one shared canvas and passes `letter="a"` to `header()` — the scripts
here do not change.  Stacked at full width they come to 171 mm, inside the
247 mm text height, so unlike Figure 2 this one does fit on a page.
