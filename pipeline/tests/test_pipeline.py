"""Unit tests for the pure-logic pipeline functions (no trajectory / file I/O).

Run from the project root with either:
    python -m unittest discover -s tests
    python -m pytest tests/            # if pytest is installed
"""
import unittest

import numpy as np

from pipeline import interface_descriptors, compute_qc, interface_metrics, exports, aggregate, \
    cdr_contacts, annotate_tcr, antigen, docking_angle, consensus


class TestAntigen(unittest.TestCase):
    def test_category(self):
        c = antigen.category
        # viral keywords win even when the organism name contains "human"
        self.assertEqual(c("Influenza A virus"), "viral")
        self.assertEqual(c("Human immunodeficiency virus 1"), "viral")
        self.assertEqual(c("human gammaherpesvirus 4"), "viral")
        self.assertEqual(c("Severe acute respiratory syndrome coronavirus 2"), "viral")
        self.assertEqual(c("Mycobacterium tuberculosis H37Rv"), "bacterial")
        self.assertEqual(c("Homo sapiens"), "tumor/self")
        self.assertEqual(c("Mus musculus"), "murine")
        self.assertEqual(c("synthetic construct"), "synthetic")
        self.assertIsNone(c(None))
        self.assertIsNone(c(""))
        # no organism, but an explicit self-peptide name -> tumor/self
        self.assertEqual(c(None, "self-peptide P1049"), "tumor/self")

    def test_extract_drops_sequence_name(self):
        # RCSB sometimes sets the peptide description to the sequence itself
        meta = {"peptide_seq": "GILGFVFTL", "rcsb": {"entities": [
            {"role": "peptide", "description": "GILGFVFTL", "organism": "Homo sapiens"}]}}
        ag = antigen.extract(meta)
        self.assertIsNone(ag["name"])
        self.assertEqual(ag["organism"], "Homo sapiens")
        self.assertEqual(ag["category"], "tumor/self")

    def test_extract_named(self):
        meta = {"peptide_seq": "LLFGYPVYV", "rcsb": {"entities": [
            {"role": "peptide", "description": "TAX PEPTIDE",
             "organism": "Human T-cell leukemia virus type I"}]}}
        ag = antigen.extract(meta)
        self.assertEqual(ag["name"], "TAX PEPTIDE")
        self.assertEqual(ag["category"], "viral")

    def test_category_name_keywords(self):
        c = antigen.category
        # RCSB leaves peptide organism empty; the antigen identity is in the name
        self.assertEqual(c(None, "TAX PEPTIDE P6A"), "viral")
        self.assertEqual(c(None, "EBV peptide LPEPLPQGQLTAY"), "viral")
        self.assertEqual(c(None, "HCMV pp65 fragment 495-503 (NLVPMVATV)"), "viral")
        self.assertEqual(c(None, "peptide from Trans-activator protein BZLF1"), "viral")
        self.assertEqual(c(None, "HPVG peptide from Epstein-Barr nuclear antigen 1"), "viral")
        self.assertEqual(c(None, "Nucleocapsid"), "viral")
        self.assertEqual(c(None, "Cancer/testis antigen 1B"), "tumor/self")  # NY-ESO-1
        self.assertEqual(c(None, "Mart-1 (26-35) peptide"), "tumor/self")
        self.assertEqual(c(None, "Melanoma peptide L7A"), "tumor/self")
        self.assertEqual(c(None, "HuD peptide"), "tumor/self")
        self.assertEqual(c(None, "EEYLKAWTF, mimotope peptide"), "synthetic")

    def test_category_name_beats_organism(self):
        # an EBV peptide annotated (wrongly) as Homo sapiens is still viral
        self.assertEqual(antigen.category("Homo sapiens", "Epstein Barr Virus peptide"), "viral")

    def test_category_host_fallback(self):
        c = antigen.category
        # no organism, no informative name -> bucket by the MHC host organism
        self.assertEqual(c(None, "NATURALLY PROCESSED OCTAPEPTIDE PBM1", "Mus musculus"), "murine")
        self.assertEqual(c(None, None, "Homo sapiens"), "tumor/self")
        self.assertIsNone(c(None, None, None))

    def test_extract_host_fallback(self):
        # mouse self-peptide: empty peptide organism, name has no pathogen/tumour
        # keyword -> falls back to the MHC entity's organism (Mus musculus -> murine)
        meta = {"peptide_seq": "INFDFNTI", "rcsb": {"entities": [
            {"role": "mhc", "description": "H-2Kb", "organism": "Mus musculus"},
            {"role": "peptide", "description": "NATURALLY PROCESSED OCTAPEPTIDE PBM1",
             "organism": None}]}}
        ag = antigen.extract(meta)
        self.assertEqual(ag["organism"], None)
        self.assertEqual(ag["category"], "murine")


class TestInterfaceDescriptors(unittest.TestCase):
    """The three qualitative recognition axes were withdrawn (their bins were terciles
    of this dataset's own run3 distribution, so a label asserted rank-in-collection,
    not a property of the molecule). What is published now is the measured ratios."""

    def _an(self, pr, hr, ac, occ=10.0, pairs=120, cdr_occ=3.0):
        return {"tcr_cdr": {"peptide_recognition_ratio": pr, "hla_restriction_ratio": hr,
                            "alpha_contribution": ac, "beta_contribution": round(1.0 - ac, 4),
                            "occupancy_total": occ, "n_tcr_phla_pairs": pairs,
                            "by_region": {"cdr1": {"peptide": {"occ": 0.0}, "mhc": {"occ": 0.0}},
                                          "cdr2": {"peptide": {"occ": 0.0}, "mhc": {"occ": 0.0}},
                                          "cdr3": {"peptide": {"occ": cdr_occ}, "mhc": {"occ": 0.0}}}}}

    def test_emits_measured_ratios(self):
        d = interface_descriptors.describe(self._an(0.4012345, 0.1, 0.7))
        self.assertTrue(d["cdr_decomposition_reliable"])
        self.assertEqual(d["peptide_recognition_ratio"], 0.4012)   # rounded, not binned
        self.assertEqual(d["hla_restriction_ratio"], 0.1)
        self.assertEqual(d["alpha_contribution"], 0.7)
        # no categorical axis survives
        for gone in ("recognition_target", "chain_bias", "interface_dynamics", "basis", "soft"):
            self.assertNotIn(gone, d)

    def test_unreliable_withholds_ratios(self):
        # occupancy exists but no CDR loop matched: the ratios are 0 by parse failure,
        # so they must be withheld rather than published as measurements. (6vma was the
        # historical example; it became reliable after the 2026-08-25 frame repair.)
        d = interface_descriptors.describe(self._an(0.0, 0.0, 0.0, occ=16.4, cdr_occ=0.0))
        self.assertFalse(d["cdr_decomposition_reliable"])
        self.assertEqual(d["cdr_unreliable_reason"], "degenerate")
        self.assertNotIn("peptide_recognition_ratio", d)

    def test_missing_ratio_is_null_not_dropped(self):
        d = interface_descriptors.describe(self._an(0.4, None, 0.7))
        self.assertIsNone(d["hla_restriction_ratio"])

    def test_none_without_tcr(self):
        self.assertIsNone(interface_descriptors.describe({"angle": {}}))


class TestCdrReliability(unittest.TestCase):
    def _br(self, cdr_occ):  # by_region with the given total CDR occupancy on peptide
        return {"cdr3": {"peptide": {"occ": cdr_occ}, "mhc": {"occ": 0.0}},
                "cdr1": {"peptide": {"occ": 0.0}, "mhc": {"occ": 0.0}},
                "cdr2": {"peptide": {"occ": 0.0}, "mhc": {"occ": 0.0}}}

    def test_reliable(self):
        self.assertIsNone(cdr_contacts.is_unreliable(
            {"occupancy_total": 10.0, "n_tcr_phla_pairs": 120, "by_region": self._br(3.0)}))

    def test_degenerate(self):
        # occupancy exists but zero CDR occupancy -> parse failure
        self.assertEqual(cdr_contacts.is_unreliable(
            {"occupancy_total": 16.4, "n_tcr_phla_pairs": 133, "by_region": self._br(0.0)}),
            "degenerate")

    def test_pair_count_outlier(self):
        # implausibly many TCR-pHLA pairs -> contact-stage artifact
        self.assertEqual(cdr_contacts.is_unreliable(
            {"occupancy_total": 50.0, "n_tcr_phla_pairs": 5105, "by_region": self._br(2.0)}),
            "pair_count_outlier")


class TestDockingPolarity(unittest.TestCase):
    # MHC platform in z=0 spread along X; peptide along +X (N->C); the inter-domain
    # vector's X-sign sets the docking polarity. Forward and reversed differ only in
    # that sign, so the crossing MAGNITUDE is identical and only polarity flips.
    _PLATFORM = [[0, 2, 0], [0, -2, 0], [30, 2, 0], [30, -2, 0]]  # symmetric -> long axis = X
    _PEP = [[5, 0, 2], [25, 0, 2]]   # pv = +X

    def _angles(self, va_pt, vb_pt):
        X = np.array(self._PLATFORM + self._PEP + [va_pt, vb_pt], dtype=float)
        return docking_angle._frame_angles(X, [0, 1, 2, 3], [6], [7], [4, 5])

    def test_forward_vs_reversed(self):
        cf, _, pf = self._angles([15, 0, 10], [25, 5, 10])   # Vβ−Vα = (+10, 5, 0)
        cr, _, prr = self._angles([15, 0, 10], [5, 5, 10])   # Vβ−Vα = (−10, 5, 0)
        self.assertEqual(pf, 1.0)               # forward polarity
        self.assertEqual(prr, -1.0)             # reversed polarity
        self.assertAlmostEqual(cf, cr, places=6)  # same crossing magnitude

    def test_near_perpendicular_undefined(self):
        # inter-domain vector almost along Y (perpendicular to the groove) -> sign noise
        _, _, p = self._angles([15, 0, 10], [15.3, 12, 10])
        self.assertEqual(p, 0.0)

    def test_directed_crossing_convention(self):
        # compute() reflects reverse-polarity frames past 90 deg (180 - acute) to build
        # the directed 0-180 crossing; forward frames keep their acute value. This locks
        # that transform (the same np.where used in docking_angle.compute).
        cf, _, pf = self._angles([15, 0, 10], [25, 5, 10])    # forward
        cr, _, prr = self._angles([15, 0, 10], [5, 5, 10])    # reversed (same magnitude)
        directed = lambda c, p: (180.0 - c) if p < 0 else c
        self.assertAlmostEqual(directed(cf, pf), cf, places=6)         # forward unchanged
        self.assertAlmostEqual(directed(cr, prr), 180.0 - cr, places=6)  # reversed -> >90
        self.assertGreater(directed(cr, prr), 90.0)


class TestReplicaStats(unittest.TestCase):
    def test_single_replica_undefined(self):
        st = aggregate._replica_stat([45.0], aggregate.REPLICA_TOL["crossing"])
        self.assertEqual(st["n"], 1)
        self.assertEqual(st["spread"], 0.0)
        self.assertIsNone(st["reproducible"])   # no second run to compare

    def test_two_replicas_within_and_over_tol(self):
        tol = aggregate.REPLICA_TOL["crossing"]  # 10 deg
        self.assertTrue(aggregate._replica_stat([45.0, 47.0], tol)["reproducible"])
        self.assertFalse(aggregate._replica_stat([30.0, 65.0], tol)["reproducible"])

    def test_three_replicas_range_based(self):
        # forward-compatible with parallel1/run1: spread/flag use the full range
        st = aggregate._replica_stat([44.0, 46.0, 52.0], aggregate.REPLICA_TOL["crossing"])
        self.assertEqual(st["n"], 3)
        self.assertEqual(st["mean"], round((44 + 46 + 52) / 3, 4))
        self.assertEqual(st["spread"], 4.0)          # (52-44)/2
        self.assertTrue(st["reproducible"])           # range 8 <= 10
        self.assertFalse(aggregate._replica_stat([44.0, 46.0, 60.0],
                                                 aggregate.REPLICA_TOL["crossing"])["reproducible"])

    def test_none_when_empty(self):
        self.assertIsNone(aggregate._replica_stat([None, None], 10.0))

    def test_ratio_tolerances_exceed_measured_replica_spread(self):
        # The reproducibility flags are only meaningful if the tolerance is wider than
        # the inter-replica scatter actually measured on the 90 complexes untouched by
        # the 2026-08-25 frame repair (medians 0.039 / 0.044 / 0.084) and over the
        # full 245-complex library (0.044 / 0.048 / 0.076, what fig4c plots).
        for key, measured_median in (("pep_recog", 0.044), ("hla_restr", 0.048),
                                     ("alpha_contrib", 0.084)):
            self.assertGreater(aggregate.REPLICA_TOL[key], measured_median)

    def test_dominant_population_not_aggregated(self):
        # The clustering descriptor was withdrawn: its inter-replica range was 29.6
        # points (34.6% of the library range), so it never entered the release.
        self.assertNotIn("dom_pop", aggregate.REPLICA_TOL)
        self.assertFalse(hasattr(aggregate, "_dyn_bin"))


class TestComputeQC(unittest.TestCase):
    def test_ok_vs_outlier(self):
        t = list(range(0, 200000, 200))
        flat = compute_qc.assess(t, [0.3] * len(t))
        self.assertEqual(flat["equilibration"], "ok")
        high = compute_qc.assess(t, [1.5] * len(t))
        self.assertEqual(high["equilibration"], "outlier")

    def test_drifting(self):
        t = list(range(0, 200000, 200))
        rising = [0.2 + 0.000004 * x for x in t]  # climbs ~0.8 nm over the run
        self.assertEqual(compute_qc.assess(t, rising)["equilibration"], "drifting")

    def test_short(self):
        self.assertIsNone(compute_qc.assess([0, 1], [0.1, 0.1]))


class TestInterfaceMetrics(unittest.TestCase):
    def test_contact_entropy(self):
        rows = [{"contact_frequency": "1.0"}, {"contact_frequency": "1.0"},
                {"contact_frequency": "1.0"}, {"contact_frequency": "1.0"}]
        e = interface_metrics.contact_entropy(rows)
        self.assertEqual(e["n_pairs"], 4)
        self.assertAlmostEqual(e["entropy_normalized"], 1.0, places=3)  # uniform -> max entropy
        # concentrated -> lower normalised entropy
        rows2 = [{"contact_frequency": "0.97"}, {"contact_frequency": "0.01"}, {"contact_frequency": "0.02"}]
        self.assertLess(interface_metrics.contact_entropy(rows2)["entropy_normalized"], 0.7)

    def test_entropy_needs_two(self):
        self.assertIsNone(interface_metrics.contact_entropy([{"contact_frequency": "1.0"}]))


class TestCDRResids(unittest.TestCase):
    def test_locate(self):
        seq = "KEVEQDRGSQSXYZAVTTDSWGKLQ"
        residues = [(i + 1, aa) for i, aa in enumerate(seq)]
        cmap = cdr_contacts.cdr_resids_for_chain(residues, {"cdr1": "DRGSQS", "cdr3": "AVTTDSWGKLQ"})
        self.assertEqual(min(cmap["cdr1"]), 6)   # 'D' is the 6th residue (1-based)
        self.assertEqual(len(cmap["cdr3"]), 11)
        self.assertEqual(cmap["_unmatched"], [])

    def test_unmatched(self):
        residues = [(i + 1, aa) for i, aa in enumerate("AAAAAA")]
        cmap = cdr_contacts.cdr_resids_for_chain(residues, {"cdr3": "WWW"})
        self.assertIn("cdr3", cmap["_unmatched"])


class TestAnnotateTCR(unittest.TestCase):
    def test_position_and_type(self):
        self.assertEqual(annotate_tcr._position("TRD"), "alpha")  # delta sits at the alpha position
        self.assertEqual(annotate_tcr._position("TRB"), "beta")
        ab = {"alpha": {"locus": "TRD"}, "beta": {"locus": "TRB"}}
        self.assertEqual(annotate_tcr.classify_type(ab), "ab")     # delta+beta = alpha/beta TCR
        gd = {"alpha": {"locus": "TRD"}, "beta": {"locus": "TRG"}}
        self.assertEqual(annotate_tcr.classify_type(gd), "gd")


class TestAggregate(unittest.TestCase):
    def test_stat(self):
        s = aggregate._stat([1, 2, 3, 4, 5])
        self.assertEqual(s["n"], 5)
        self.assertEqual(s["median"], 3)
        self.assertEqual(s["min"], 1)
        self.assertEqual(s["max"], 5)
        self.assertEqual(aggregate._stat([])["n"], 0)

    def test_stat_interpolates(self):
        # even-length input: median/quartiles must LINEARLY interpolate (numpy),
        # not fall back to the old biased nearest-rank s[int(p*n)].
        s = aggregate._stat([1, 2, 3, 4])
        self.assertEqual(s["median"], 2.5)
        self.assertEqual(s["p25"], 1.75)
        self.assertEqual(s["p75"], 3.25)

    def test_hla_group(self):
        # Grouping coarsens to the HLA locus so finely- and coarsely-annotated
        # alleles of the same locus land in one consistent bucket.
        self.assertEqual(aggregate._hla_group("HLA-A*02:01"), "HLA-A")
        self.assertEqual(aggregate._hla_group("HLA-A"), "HLA-A")
        self.assertEqual(aggregate._hla_group("HLA-E"), "HLA-E")
        self.assertIsNone(aggregate._hla_group(None))


class TestReplicaConsensus(unittest.TestCase):
    # majority_label survives only for genuinely categorical, non-derived attributes
    # (e.g. TCR chain type); the three tercile-binned recognition axes it used to vote
    # on were withdrawn in favour of numeric_consensus over the measured ratios.
    def test_majority_beats_primary(self):
        # primary (run3) is the minority label -> the majority wins
        label, agree, n = consensus.majority_label(
            ["alphabeta", "gammadelta", "gammadelta"],
            primary="alphabeta", tie_break="alphabeta")
        self.assertEqual(label, "gammadelta")
        self.assertFalse(agree)
        self.assertEqual(n, 3)

    def test_tie_breaks_to_declared_fallback(self):
        label, agree, _ = consensus.majority_label(
            ["alphabeta", "gammadelta"], primary="gammadelta", tie_break="alphabeta")
        self.assertEqual(label, "alphabeta")
        self.assertFalse(agree)

    def test_unanimous_agrees(self):
        label, agree, n = consensus.majority_label(["alphabeta", "alphabeta"])
        self.assertEqual(label, "alphabeta")
        self.assertTrue(agree)

    def test_single_replica_agreement_none(self):
        label, agree, n = consensus.majority_label(["alphabeta"])
        self.assertEqual(label, "alphabeta")
        self.assertIsNone(agree)
        self.assertEqual(n, 1)

    def test_numeric_consensus_reports_spread(self):
        st = consensus.numeric_consensus([0.40, 0.44, 0.42])
        self.assertEqual(st["n"], 3)
        self.assertEqual(st["mean"], 0.42)
        self.assertEqual(st["min"], 0.40)
        self.assertEqual(st["max"], 0.44)
        self.assertEqual(st["range"], 0.04)
        self.assertAlmostEqual(st["sd"], 0.02, places=4)

    def test_numeric_consensus_single_and_empty(self):
        st = consensus.numeric_consensus([0.31, None])
        self.assertEqual((st["n"], st["mean"], st["sd"], st["range"]), (1, 0.31, 0.0, 0.0))
        self.assertIsNone(consensus.numeric_consensus([None, None]))

    def test_all_none_returns_none(self):
        self.assertEqual(consensus.majority_label([None, None]), (None, None, 0))

    def test_bool_majority_and_tie_weight(self):
        # strict majority of booleans
        v, agree, _ = consensus.majority_bool([True, True, False])
        self.assertTrue(v)
        self.assertFalse(agree)
        # even split -> decided by mean reversed fraction
        v2, _, _ = consensus.majority_bool([True, False], weights=[0.9, 0.8])
        self.assertTrue(v2)
        v3, _, _ = consensus.majority_bool([True, False], weights=[0.1, 0.2])
        self.assertFalse(v3)

    def test_qc_worst_case(self):
        # a single outlier replica makes the complex non-ok (never hidden)
        r = consensus.qc_consensus(["ok", "ok", "outlier"])
        self.assertEqual(r["consensus"], "outlier")
        self.assertFalse(r["reproducible"])
        self.assertEqual(r["counts"], {"ok": 2, "drifting": 0, "outlier": 1})

    def test_qc_unanimous_and_single(self):
        self.assertTrue(consensus.qc_consensus(["ok", "ok"])["reproducible"])
        self.assertIsNone(consensus.qc_consensus(["ok"])["reproducible"])
        self.assertIsNone(consensus.qc_consensus([None])["consensus"])


class TestExports(unittest.TestCase):
    def test_feature_table(self):
        meta = {"traj_id": "x_run3", "pdb_id": "x", "peptide_seq": "SIINFEKL", "peptide_length": 8,
                "tcr": {"chains": {"alpha": {"v_gene": "TRAV1", "cdr3": "CAVR"}}}}
        analysis = {"tcr_cdr": {"peptide_recognition_ratio": 0.3},
                    "geometry": {"bulge_height_angstrom": 5.0}}
        csv = exports.feature_table_csv(meta, analysis)
        self.assertIn("category,feature,value", csv)
        self.assertIn("SIINFEKL", csv)
        self.assertIn("peptide_recognition_ratio,0.3", csv)
        self.assertIn("TRAV1", csv)


if __name__ == "__main__":
    unittest.main()
