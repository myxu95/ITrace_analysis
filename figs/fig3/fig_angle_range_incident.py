#!/usr/bin/env python
"""Incident angle: a median line over each complex's fluctuation band, all 245.

The incident angle is the tilt at which the TCR approaches the pMHC surface.
Where the crossing angle says how the receptor is rotated over the groove, the
incident angle says how steeply it sits on it, so the two are independent
descriptions of the same approach and each gets its own figure.

This is a system-level statistic, not a validation panel: the band is what one
complex's angle does over 600 ns of sampling, and reading the band across the
rank axis gives the range of approach geometry the library covers.
The three replicas are pooled into one 3 x 1001-frame sample per complex,
because here they are three views of the same starting structure rather than
three things to compare.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from anglestrip import PH, strip
from figstyle import *
from repdata import N_TOTAL, library, quant

LEGEND = {}

VMAX = 70.0
TICKS = [0, 10, 20, 30, 40, 50, 60, 70]


def draw(ax, top, w):
    rows = library()
    med = [r["inc_p50"] for r in rows]
    span = [r["inc_p95"] - r["inc_p5"] for r in rows]
    narrow = min(rows, key=lambda r: r["inc_p95"] - r["inc_p5"])
    wide = max(rows, key=lambda r: r["inc_p95"] - r["inc_p5"])
    n_wide = sum(1 for s in span if s > 20)

    LEGEND.update(
        title="Incident angle of the TCR approach",
        subtitle=f"All {N_TOTAL} complexes, ranked by median. The incident angle is the "
                 f"tilt of the TCR variable domains against the pMHC surface. Line, the "
                 f"median of one complex over its pooled 3 × 1001 frames; shading, that "
                 f"complex\u2019s 5th to 95th percentile, so the height of the band is how far "
                 f"one complex\u2019s angle fluctuates and the rise of the line is the range "
                 f"the library covers. Percentiles rather than the full range, because over "
                 f"three thousand frames the extremes describe a single frame.",
        notes=[f"Median incident angle runs from {min(med):.1f}° to {max(med):.1f}° across "
               f"the library, with the middle half of complexes between "
               f"{quant(med, .25):.1f}° and {quant(med, .75):.1f}°.",
               f"A complex occupies a {quant(span, .5):.1f}° interval at the median "
               f"({quant(span, .25):.1f}–{quant(span, .75):.1f}° over the middle half), "
               f"narrow against the {max(med) - min(med):.0f}° the library spans, so for "
               f"most complexes the approach tilt is set by the complex rather than "
               f"sampled over during the run. That does not hold in the tail: {n_wide} "
               f"complexes exceed 20°, where the band records genuine displacement of the "
               f"TCR over the groove — progressive within one replica, or a persistent "
               f"offset between replicas.",
               f"The narrowest interval is {narrow['pdb_id']} at "
               f"{narrow['inc_p95'] - narrow['inc_p5']:.1f}° and the widest "
               f"{wide['pdb_id']} at {wide['inc_p95'] - wide['inc_p5']:.1f}°."],
    )
    return strip(ax, top, w, rows, "inc", VMAX, TICKS,
                 "Incident angle (°)", LEGEND)


if __name__ == "__main__":
    render(draw, "fig_angle_range_incident")
