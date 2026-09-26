"""Per-trajectory interface descriptors — quantitative values only.

Publishes the measured CDR-decomposition ratios directly, instead of binning them
into dataset-relative categorical labels.

**Why this module no longer classifies (2026-08-25).** It previously emitted a
three-axis "dynamic recognition profile" — ``recognition_target`` (peptide-centered /
intermediate / HLA-centered), ``chain_bias`` (alpha-dominant / balanced /
beta-dominant) and ``interface_dynamics`` (rigid-lock / intermediate / dynamic) — by
cutting each scalar at terciles of this dataset's own run3 distribution. A label such
as "peptide-centered" therefore asserted only *"in the lower third of these 245
systems"*: a statement about a complex's rank inside our collection, not a property of
the molecule, and one that would change if the collection changed. Those bins are
withdrawn; the scalars they were derived from are published instead, so downstream
users can apply whatever threshold their own question requires.

``interface_dynamics`` is withdrawn along with its underlying scalar
(``dominant_population_percent``, the dominant interface-cluster population). That
quantity is not reproducible across independent replicas — median inter-replica range
29.6 percentage points, 34.6% of the library range, measured on 90 complexes whose
three replicas were all untouched by the 2026-08-25 frame repair — because the
upstream clustering cuts a dendrogram at a fixed absolute cutoff over a
max-normalised distance matrix, making the cluster count hostage to a single outlier
frame pair. No interface-clustering descriptor is served.

Published per trajectory:
  * ``peptide_recognition_ratio``  — CDR3 -> peptide share of TCR interface contacts
  * ``hla_restriction_ratio``      — CDR1/2 -> HLA share
  * ``alpha_contribution``         — alpha-chain share of TCR interface contacts
  * ``cdr_decomposition_reliable`` — False when the decomposition is degenerate or
    artefactual; the three ratios are then withheld and the reason is given.

These three ratios ARE reproducible, in the sense the withdrawn one was not. Over the
full 245-complex library their median inter-replica ranges are 0.044 / 0.048 / 0.076 --
10-12% of the library range for each -- and only 2 complexes (both on
``hla_restriction_ratio``) span more than half the library range, against 22 of 90 for
the withdrawn dominant-population descriptor. On the same 90 untouched complexes the
medians are 0.039 / 0.044 / 0.084, so the repair is not what produces this.

They are not uniformly tight, and the served data says so rather than hiding it: 24 /
39 / 89 complexes of 245 exceed the 0.10 tolerance (``frac_reproducible`` 0.90 / 0.84 /
0.64). ``alpha_contribution`` is the softest of the three. Every served value carries
its measured spread, so a user can apply their own bar.

Complex-level aggregation is a mean plus the *measured* inter-replica spread
(``consensus.numeric_consensus``), which replaces the old boolean ``*_agree`` flags:
a reported spread tells a user how far the three independent runs actually
disagreed, where the flag only said whether they landed in the same bin.

    IMMUNO_WEB_DATA=.../immuno-dyn python -m pipeline.interface_descriptors
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from . import config
from . import cdr_contacts

#: The quantitative descriptors served for every trajectory with a reliable
#: CDR decomposition. Order is the served/report order.
RATIO_KEYS = ("peptide_recognition_ratio", "hla_restriction_ratio", "alpha_contribution")

#: Key written into meta.json. (Was ``recognition_mode`` before 2026-08-25.)
META_KEY = "interface_descriptors"

#: Library-level summary file. (Was ``recognition_modes.json``.)
SUMMARY_NAME = "interface_descriptors.json"


def describe(analysis: dict) -> dict | None:
    """Quantitative interface descriptors for one trajectory's ``analysis.json``.

    Returns ``None`` when the trajectory has no CDR decomposition at all. Otherwise
    returns the reliability flag, plus the three ratios when they can be trusted.
    Unreliable decompositions keep the flag and carry the reason, so a consumer can
    tell "not measurable here" apart from "not measured".
    """
    tcr = analysis.get("tcr_cdr")
    if not tcr:
        return None

    reason = cdr_contacts.is_unreliable(tcr)
    out: dict = {"cdr_decomposition_reliable": reason is None}
    if reason is not None:
        out["cdr_unreliable_reason"] = reason
        return out

    for key in RATIO_KEYS:
        value = tcr.get(key)
        out[key] = round(float(value), 4) if value is not None else None
    return out


def _summarize(values: list[float]) -> dict:
    """Distribution summary for one descriptor across the library."""
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return {"n": 0}
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0

    def q(p: float) -> float:
        pos = p * (n - 1)
        lo = math.floor(pos)
        hi = math.ceil(pos)
        return vals[lo] if lo == hi else vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)

    return {
        "n": n,
        "mean": round(mean, 4),
        "sd": round(math.sqrt(var), 4),
        "min": round(vals[0], 4),
        "q1": round(q(0.25), 4),
        "median": round(q(0.50), 4),
        "q3": round(q(0.75), 4),
        "max": round(vals[-1], 4),
    }


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", nargs="*", help="limit to these trajectory ids")
    args = ap.parse_args(argv)

    root = Path(config.WEB_DATA)
    written = unreliable = skipped = 0
    collected: dict[str, list[float]] = {k: [] for k in RATIO_KEYS}

    for meta_path in sorted(root.glob("*/meta.json")):
        traj_id = meta_path.parent.name
        if args.ids and traj_id not in args.ids:
            continue
        an_path = meta_path.parent / "analysis" / "analysis.json"
        if not an_path.exists():
            skipped += 1
            continue
        analysis = json.loads(an_path.read_text())
        desc = describe(analysis)
        meta = json.loads(meta_path.read_text())
        # Drop the withdrawn categorical block wherever it survives from an
        # older build, so a partial re-run cannot leave stale labels on disk.
        meta.pop("recognition_mode", None)
        if desc is None:
            meta.pop(META_KEY, None)
            skipped += 1
        else:
            meta[META_KEY] = desc
            written += 1
            if desc["cdr_decomposition_reliable"]:
                for key in RATIO_KEYS:
                    if desc.get(key) is not None:
                        collected[key].append(desc[key])
            else:
                unreliable += 1
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    summary = {
        "n_trajectories": written,
        "n_unreliable_cdr": unreliable,
        "descriptors": {k: _summarize(v) for k, v in collected.items()},
        "note": ("Quantitative interface descriptors. No categorical recognition-mode "
                 "labels are derived from these values; no interface-clustering "
                 "descriptor is served."),
    }
    out_dir = root
    (out_dir / SUMMARY_NAME).write_text(json.dumps(summary, indent=2) + "\n")
    # Remove the superseded label-distribution summary if a previous build left it.
    legacy = out_dir / "recognition_modes.json"
    if legacy.exists():
        legacy.unlink()

    print(f"written={written} unreliable_cdr={unreliable} skipped={skipped}")
    for key, stats in summary["descriptors"].items():
        print(f"  {key:28s} {stats}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
