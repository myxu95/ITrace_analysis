#!/usr/bin/env python
"""Crossing angle: a median line over each complex's fluctuation band, all 245.

The crossing angle is the rotation of the TCR variable domains over the axis of
the peptide-binding groove, measured as a directed angle on 0-180 so that the
two docking polarities stay distinguishable instead of being folded onto each
other by a magnitude.

A system-level statistic: the band is one complex's own sampling, and the rank
axis gives the range of docking rotation the library covers.  Replicas are
pooled, since here they are three views of one starting structure.

Four complexes dock in reverse polarity and sit above 120 degrees.  They are
drawn in the contrast colour rather than dropped: the polarity is real geometry
reported by the pipeline for every trajectory of those four, not an artefact,
and an atlas that hides its four unusual entries is not an atlas.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from anglestrip import PH, strip
from figstyle import *
from repdata import N_TOTAL, library, quant

LEGEND = {}

VMAX = 180.0
TICKS = [0, 30, 60, 90, 120, 150, 180]
REVERSED = ("5sws", "5swz", "7jwi", "9gv7")     # angle.reversed_polarity, all frames


def draw(ax, top, w):
    rows = library()
    fwd = [r for r in rows if r["pdb_id"] not in REVERSED]
    med = [r["cross_p50"] for r in fwd]
    span = [r["cross_p95"] - r["cross_p5"] for r in fwd]

    LEGEND.update(
        title="Crossing angle of the TCR over the groove",
        subtitle=f"All {N_TOTAL} complexes, ranked by median. The crossing angle is "
                 f"directed on 0–180°, so the {len(REVERSED)} complexes that dock in "
                 f"reverse polarity ({', '.join(p.upper() for p in REVERSED)}) stay separable, and they are "
                 f"drawn as their own band in the contrast colour rather than joined to "
                 f"the rest. Line, the median of one complex over its pooled "
                 f"3 × 1001 frames; shading, that complex\u2019s 5th to 95th percentile.",
        notes=[f"Across the {len(fwd)} forward-polarity complexes the median crossing "
               f"angle runs from {min(med):.0f}° to {max(med):.0f}°, with the middle half "
               f"between {quant(med, .25):.0f}° and {quant(med, .75):.0f}°.",
               f"A complex occupies a {quant(span, .5):.0f}° interval at the median "
               f"({quant(span, .25):.0f}–{quant(span, .75):.0f}° over the middle half), "
               f"against the {max(med) - min(med):.0f}° the forward-polarity library "
               f"spans, so docking rotation separates complexes about as cleanly as "
               f"approach tilt does.",
               f"The {len(REVERSED)} reverse-polarity complexes carry reversed_polarity "
               f"in every one of their {len(REVERSED) * 3} trajectories, so their position "
               f"above 120° is the geometry of the deposited complex and not a "
               f"trajectory that turned over."],
    )
    return strip(ax, top, w, rows, "cross", VMAX, TICKS,
                 "Crossing angle (°)", LEGEND, exception=REVERSED,
                 classes=("Forward polarity", "Reverse polarity"))


if __name__ == "__main__":
    render(draw, "fig_angle_range_crossing")
