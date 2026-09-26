#!/usr/bin/env python
"""Figure 3 -- what the trajectories say about the peptide and the docking, as
one composed plate.

Four panels in two rows, which is the arrangement the manuscript already
carries.  The top row is the peptide: how far it bulges out of the groove
against its length, and where along it the TCR actually sits.  The bottom row
is the docking geometry, the same drawing of two angles -- ranked band, median
line -- so the second panel is read for free once the first is.

c and d already share their whole drawing (anglestrip.py), so their frames
align by construction.  a and b cannot: b is two stacked box rows and a is one,
and forcing a common plot height would only stretch a to no purpose.  They hang
from the same top edge, which is the alignment that means something here.

Nothing is redrawn.  Each panel is the code that renders it on its own, called
into a sub-axes whose coordinates are still millimetres, so the plate is vector
throughout.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from compose import compose

import fig_angle_range_crossing as crossing
import fig_angle_range_incident as incident
import fig_bulge_by_length as bulge
import fig_peptide_position as position

# b carries nine positions in two rows and needs the width more than a, whose
# axis is seven lengths with an elision in it.
ROWS = [[("a", bulge.draw, 0.94), ("b", position.draw, 1.06)],
        [("c", crossing.draw), ("d", incident.draw)]]

if __name__ == "__main__":
    compose("figure3", ROWS, gx=6.0, gy=5.0,
            out=Path(__file__).resolve().parent)
