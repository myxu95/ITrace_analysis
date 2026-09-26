"""Unit tests for the concerted-motion module's pure logic (no trajectory I/O).

The DCCM / network compute path needs mdtraj + a real trajectory; these tests
cover the pure helpers and the math contract instead, and skip cleanly if the
optional deps (networkx/matplotlib) are not installed.

Run from the project root:
    python -m unittest discover -s tests
"""
import unittest

import numpy as np

try:
    import networkx as nx
    from pipeline import concerted_motion as cm
    _HAVE = True
except Exception:  # missing networkx/matplotlib -> skip rather than break discovery
    _HAVE = False


@unittest.skipUnless(_HAVE, "concerted_motion deps (networkx/matplotlib) not installed")
class TestConcertedMotion(unittest.TestCase):
    def test_mhc_chain_is_largest_remaining(self):
        meta = {"chains": [{"id": "A", "n_residues": 275}, {"id": "B", "n_residues": 100},
                           {"id": "C", "n_residues": 9}, {"id": "D", "n_residues": 115},
                           {"id": "E", "n_residues": 209}]}
        # peptide=C, TCRa=D, TCRb=E -> MHC is the largest of {A, B} = A
        self.assertEqual(cm._mhc_chain(meta, "C", "D", "E"), "A")

    def test_mhc_chain_none_when_no_candidate(self):
        meta = {"chains": [{"id": "C", "n_residues": 9}]}
        self.assertIsNone(cm._mhc_chain(meta, "C", "D", "E"))

    def test_resolve_cdr3_exact_subsequence(self):
        # chain D residues spell 'KAVTTDSWGKLQ'; CDR3a 'AVTTDSWGKLQ' starts at offset 1
        resn = ["LYS", "ALA", "VAL", "THR", "THR", "ASP", "SER", "TRP", "GLY", "LYS", "LEU", "GLN"]
        chain = ["D"] * len(resn)
        nodes = cm._resolve_cdr3_nodes(chain, resn, "D", "E", "AVTTDSWGKLQ", None)
        self.assertEqual(nodes, list(range(1, 12)))

    def test_resolve_cdr3_absent_returns_empty(self):
        resn = ["GLY", "GLY", "ALA"]
        nodes = cm._resolve_cdr3_nodes(resn and ["D"] * 3, resn, "D", "E", "WWWW", None)
        self.assertEqual(nodes, [])

    def test_dccm_symmetric_unit_diag_in_range(self):
        # replicate the module's DCCM math and assert its invariants
        rng = np.random.default_rng(0)
        disp = rng.standard_normal((50, 12, 3))
        disp -= disp.mean(0, keepdims=True)
        F = disp.shape[0]
        cov = np.einsum("fix,fjx->ij", disp, disp) / F
        msf = np.diag(cov).copy()
        den = np.sqrt(np.outer(msf, msf)); den[den == 0] = 1.0
        dccm = np.clip(cov / den, -1.0, 1.0)
        self.assertTrue(np.allclose(dccm, dccm.T, atol=1e-12))
        self.assertTrue(np.allclose(np.diag(dccm), 1.0, atol=1e-9))
        self.assertGreaterEqual(dccm.min(), -1.0)
        self.assertLessEqual(dccm.max(), 1.0)

    def test_networkx_3_6_api_present(self):
        # env-regression guard for the three networkx calls compute() relies on
        G = nx.path_graph(6)
        for u, v in G.edges():
            G[u][v]["weight"] = 1.0
            G[u][v]["length"] = 0.5
        comms = nx.community.louvain_communities(G, weight="weight", seed=42)
        self.assertIsInstance(nx.community.modularity(G, comms, weight="weight"), float)
        self.assertIsInstance(nx.betweenness_centrality(G, weight="length", normalized=True), dict)

    def test_analysis_key(self):
        self.assertEqual(cm.ANALYSIS_KEY, "concerted_motion")
        self.assertEqual(cm.SIDECAR, "concerted_motion.json")

    def test_ranges_compress(self):
        # contiguous runs collapse; gaps split; unsorted/dupes handled
        self.assertEqual(cm._ranges([3, 1, 2, 5]), [[1, 3], [5, 5]])
        self.assertEqual(cm._ranges([10, 11, 11, 12]), [[10, 12]])

    def test_community_blocks(self):
        # comm 0: chain A resSeq 5,6,8 (gap at 7) + 2 MHC/1 PEP; comm 1: chain B 10,11
        comm_label = np.array([0, 0, 0, 1, 1])
        node_chain = ["A", "A", "A", "B", "B"]
        node_resseq = [5, 6, 8, 10, 11]
        node_role = ["MHC", "MHC", "PEP", "TCRa", "TCRa"]
        blocks = cm._community_blocks(comm_label, node_chain, node_resseq, node_role)
        self.assertEqual([b["id"] for b in blocks], [0, 1])     # size-ranked order
        self.assertEqual(blocks[0]["n"], 3)
        self.assertEqual(blocks[0]["ranges"]["A"], [[5, 6], [8, 8]])  # gap at 7 splits
        self.assertEqual(blocks[0]["roles"], {"MHC": 2, "PEP": 1})
        self.assertEqual(blocks[1]["ranges"]["B"], [[10, 11]])
        # residue count is conserved (sum of block sizes == assigned nodes)
        self.assertEqual(sum(b["n"] for b in blocks), 5)

    def test_community_blocks_skips_unassigned(self):
        blocks = cm._community_blocks(np.array([-1, 0, 0]), ["A", "A", "A"], [1, 2, 3],
                                      ["MHC", "MHC", "MHC"])
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0]["n"], 2)


if __name__ == "__main__":
    unittest.main()
