"""Clustering view - conformational states."""

from pathlib import Path
import pandas as pd

from ..base import Presenter, RenderContext
from ..loaders import read_json, read_csv
from ..formatters import pct, num
from .. import register_view


@register_view("clustering")
class ClusteringPresenter(Presenter):
    """
    Clustering view - conformational states.

    Shows interface clustering results: number of clusters, dominant
    cluster, and transitions between states.
    """

    view_name = "clustering"
    max_length = 3000
    # D-B1: TICA / FEL / interface-cluster populations and assignments —
    # Complex-level dynamic (ensemble) information.
    spatial_layer = "complex"
    sub_flavors_served = ("dynamic",)

    def render(self, context: RenderContext) -> str:
        """Render clustering view."""
        locator = context.locator
        case_id = locator.get_case_id()

        lines = [f"# Clustering: {case_id}\n"]

        # Load clustering data
        cluster_root = locator.get_module_root("inter_cluster")
        if not cluster_root:
            lines.append("_Clustering analysis not found. Run clustering analysis first._")
            return "\n".join(lines)

        summary_path = cluster_root / "analysis/conformation/interface_clustering/interface_clustering_summary.json"
        summary = read_json(summary_path)

        if not summary:
            lines.append("_Clustering data not available._")
            return "\n".join(lines)

        # Render sections
        lines.append(self._render_summary(summary))
        lines.append(self._render_clusters(summary))
        lines.append(self._render_transitions(summary))

        # Try to add feature digest
        digest = self._render_feature_digest(cluster_root)
        if digest:
            lines.append(digest)

        # Sources
        lines.append(self._format_sources(
            summary_path.relative_to(locator.case_dir)
        ))

        output = "\n".join(lines)
        return self._cap_length(output)

    def _render_summary(self, summary: dict) -> str:
        """Render clustering summary."""
        lines = ["## Summary\n"]

        n_clusters = summary.get("n_clusters", 0)
        lines.append(f"- **Number of clusters**: {n_clusters}")

        dominant_frac = summary.get("dominant_cluster_fraction", 0)
        dominant_id = summary.get("dominant_cluster_id", 0)
        lines.append(f"- **Dominant cluster**: {dominant_id} ({pct(dominant_frac)} of frames)")

        # Convergence if available
        is_converged = summary.get("is_converged")
        if is_converged is not None:
            icon = "✅" if is_converged else "⚠️"
            lines.append(f"- **Convergence**: {icon} {'Yes' if is_converged else 'No'}")

        return "\n".join(lines) + "\n"

    def _render_clusters(self, summary: dict) -> str:
        """Render cluster details."""
        lines = ["## Clusters\n"]

        clusters = summary.get("clusters", [])
        if not clusters:
            lines.append("_Cluster details not available_")
            return "\n".join(lines)

        for cluster in clusters:
            cluster_id = cluster.get("cluster_id", "?")
            size = cluster.get("size", 0)
            fraction = cluster.get("fraction", 0)

            lines.append(f"**Cluster {cluster_id}**: {size} frames ({pct(fraction)})")

            # Add representative frame if available
            rep_frame = cluster.get("representative_frame")
            if rep_frame is not None:
                lines.append(f"  - Representative frame: {rep_frame}")

            # Add centroid info if available
            centroid = cluster.get("centroid")
            if centroid:
                lines.append(f"  - Centroid: {centroid}")

            lines.append("")

        return "\n".join(lines)

    def _render_transitions(self, summary: dict) -> str:
        """Render cluster transitions."""
        lines = ["## Transitions\n"]

        transitions = summary.get("transitions", [])
        if not transitions:
            lines.append("_Transition data not available_")
            return "\n".join(lines)

        lines.append("**Major transitions**:")
        for trans in transitions[:5]:  # Top 5 transitions
            from_cluster = trans.get("from", "?")
            to_cluster = trans.get("to", "?")
            count = trans.get("count", 0)

            lines.append(f"- Cluster {from_cluster} → {to_cluster}: {count} transitions")

        return "\n".join(lines) + "\n"

    def _render_feature_digest(self, cluster_root: Path) -> str:
        """Render feature digest if available."""
        digest_path = cluster_root / "analysis/conformation/interface_clustering/cluster_feature_digest.csv"
        df = read_csv(digest_path)

        if df.empty:
            return ""

        lines = ["## Feature digest\n"]

        # Show top features that distinguish clusters
        if "feature" in df.columns and "importance" in df.columns:
            df = df.sort_values("importance", ascending=False).head(5)

            lines.append("**Top distinguishing features**:")
            for _, row in df.iterrows():
                feature = row.get("feature", "")
                importance = row.get("importance", 0)
                lines.append(f"- {feature}: {num(importance, 3)}")

        return "\n".join(lines) + "\n"
