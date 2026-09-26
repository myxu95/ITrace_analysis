# Figure 4 — technical validation

The figure Scientific Data's mandatory **Technical Validation** section reads
on. One script builds the whole plate:

    make_fig4.py     141 x 127 mm  ->  pdf/ svg/ png/  (fig4, fig4a..fig4d)

Four panels, lettered in the order the manuscript first cites them:

| | Panel | Unit |
|---|---|---|
| a | Q, the fraction of native contacts retained | complex (245) |
| b | Incident angle of the TCR approach | complex (245) |
| c | Essential-subspace RMSIP between replicas | complex (245) |
| d | Dihedral-PCA free-energy surfaces of one released entry | one entry |

a–c are the same ranked-strip construct, so they are contiguous; d is drawn in
another register and therefore goes last. Each panel also renders on its own
(`fig4a.pdf` … `fig4d.pdf`) in case the plate is ever split.

Data comes from `../repdata.py` (`replicas`, `icc`, `quant`) and
`../landdata.py` (`exemplar`). `icc()` is one-way random-effects ICC(1,1); the
manuscript quotes the same numbers.

## Layout

Output is filed by format:

    pdf/   submission masters (vector)
    svg/   editable source
    png/   600 dpi, for looking at

That split is switched on by the `.split-by-type` marker file in this folder,
which `figstyle._dest()` looks for. It is a property of the directory, so any
script rendered here — including ones written later — files itself correctly
without being told to. **No other figure folder opts in**; delete the marker to
go back to flat output.

## Not part of the plate

`fig_interface_pca.py` and `fig_rmsd_replicas.py` are exploratory probes that
nothing imports; they are kept because the data work in them was expensive.

The v1/v2 single-question panels this figure grew out of
(`fig_replica_incident*`, `fig_replica_fnat*`, `fig_drift*`, plus their shared
`repstrip.py` / `pairscatter.py` and the `preview_v2_plate.py` sheet) now live
in `../_retired/fig4_panels/`. `_v2_archive/` is the dated review snapshot they
came from; everything in it is a copy and nothing depends on it.
