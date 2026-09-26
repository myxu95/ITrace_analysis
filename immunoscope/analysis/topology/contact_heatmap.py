"""Contact heatmap construction and plotting."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


@dataclass
class ContactHeatmapArtifacts:
    """Artifact paths for a single contact heatmap."""

    matrix_file: str
    image_file: str
    n_rows: int
    n_cols: int


class ContactHeatmapPlotter:
    """Generate contact-frequency heatmaps from contact report result tables."""

    def __init__(self, style: str = "seaborn-v0_8-white"):
        self.style = style
        self.cmap = LinearSegmentedColormap.from_list(
            "immunoscope_contact_sage",
            ["#fffaf2", "#eef3ed", "#c8d8cf", "#7da48f", "#2f5d50"],
        )

    def build_matrix(
        self,
        contacts: pd.DataFrame,
        interaction_class: str,
        tcr_chain: str | None = None,
    ) -> pd.DataFrame:
        """Build a residue-pair contact-frequency matrix."""
        if interaction_class == "groove_tcr":
            filtered = contacts[
                (contacts["interaction_class"] == "hla_tcr")
                & (contacts["mhc_region"] == "groove")
            ].copy()
        else:
            filtered = contacts[contacts["interaction_class"] == interaction_class].copy()
        if tcr_chain is not None:
            filtered = filtered[filtered["tcr_chain"] == tcr_chain].copy()
        if filtered.empty:
            return pd.DataFrame()

        filtered["x_order"] = filtered["phla_resid"].astype(int)
        filtered["y_order"] = filtered["tcr_resid"].astype(int)
        filtered["x_label"] = filtered["phla_residue"]
        filtered["y_label"] = filtered["tcr_residue"]

        pivot = filtered.pivot_table(
            index="y_label",
            columns="x_label",
            values="contact_frequency",
            aggfunc="max",
            fill_value=0.0,
        )

        x_order = (
            filtered[["x_label", "x_order"]]
            .drop_duplicates()
            .sort_values(["x_order", "x_label"])["x_label"]
            .tolist()
        )
        y_order = (
            filtered[["y_label", "y_order"]]
            .drop_duplicates()
            .sort_values(["y_order", "y_label"])["y_label"]
            .tolist()
        )

        pivot = pivot.reindex(index=y_order, columns=x_order, fill_value=0.0)
        return pivot

    def save_heatmap(
        self,
        matrix: pd.DataFrame,
        output_prefix: Path,
        title: str,
        xlabel: str,
        ylabel: str,
    ) -> ContactHeatmapArtifacts | None:
        """Save the matrix CSV and PNG heatmap."""
        if matrix.empty:
            return None

        output_prefix.parent.mkdir(parents=True, exist_ok=True)

        matrix_file = output_prefix.with_suffix(".csv")
        image_file = output_prefix.with_suffix(".png")
        matrix.to_csv(matrix_file)

        n_cols = max(len(matrix.columns), 1)
        n_rows = max(len(matrix.index), 1)
        figsize = (max(8, n_cols * 0.45), max(6, n_rows * 0.4))

        plt.style.use(self.style)
        fig, ax = plt.subplots(figsize=figsize)
        display_matrix = matrix.copy()
        display_matrix = display_matrix.where(display_matrix > 0, np.nan)

        sns.heatmap(
            display_matrix,
            ax=ax,
            cmap=self.cmap,
            vmin=0.0,
            vmax=1.0,
            cbar=True,
            cbar_kws={
                "label": "Contact Frequency",
                "ticks": [0.0, 0.25, 0.5, 0.75, 1.0],
                "shrink": 0.92,
                "pad": 0.02,
            },
            linewidths=0.25,
            linecolor="#f3f4f6",
            square=False,
        )
        ax.set_title(title, fontsize=13, pad=12)
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.tick_params(axis="x", rotation=60, labelsize=8, length=0)
        ax.tick_params(axis="y", rotation=0, labelsize=8, length=0)
        fig.tight_layout()
        fig.savefig(image_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

        return ContactHeatmapArtifacts(
            matrix_file=str(matrix_file),
            image_file=str(image_file),
            n_rows=len(matrix.index),
            n_cols=len(matrix.columns),
        )

    def _save_variant(
        self,
        contacts: pd.DataFrame,
        output_root: Path,
        interaction_class: str,
        title: str,
        xlabel: str,
        tcr_chain: str | None,
        suffix: str,
    ) -> dict:
        matrix = self.build_matrix(
            contacts=contacts,
            interaction_class=interaction_class,
            tcr_chain=tcr_chain,
        )
        artifacts = self.save_heatmap(
            matrix=matrix,
            output_prefix=output_root / self._build_output_stem(interaction_class, suffix),
            title=title,
            xlabel=xlabel,
            ylabel="TCR Residues",
        )
        if artifacts is None:
            return {
                "matrix_file": None,
                "image_file": None,
                "n_rows": 0,
                "n_cols": 0,
            }
        return {
            "matrix_file": artifacts.matrix_file,
            "image_file": artifacts.image_file,
            "n_rows": artifacts.n_rows,
            "n_cols": artifacts.n_cols,
        }


    @staticmethod
    def _build_output_stem(interaction_class: str, suffix: str) -> str:
        interaction_map = {
            "peptide_tcr": "pep",
            "hla_tcr": "hla",
            "groove_tcr": "groove",
        }
        suffix_map = {
            "alpha": "TCRa",
            "beta": "TCRb",
        }
        left = interaction_map[interaction_class]
        right = suffix_map[suffix]
        return f"{left}_{right}_heatmap"

    def generate_bundle(self, contacts: pd.DataFrame, output_dir: str | Path) -> dict:
        """Generate peptide/HLA/groove heatmaps split by TCR alpha/beta."""
        output_root = Path(output_dir)
        outputs: dict[str, dict] = {}

        specs: Iterable[tuple[str, str, str]] = (
            ("peptide_tcr", "Peptide-TCR Contact Frequency Heatmap", "Peptide Residues"),
            ("hla_tcr", "HLA-TCR Contact Frequency Heatmap", "HLA Residues"),
            ("groove_tcr", "Groove-TCR Contact Frequency Heatmap", "Groove Residues"),
        )

        for interaction_class, title, xlabel in specs:
            outputs[interaction_class] = {
                "alpha": self._save_variant(
                    contacts=contacts,
                    output_root=output_root,
                    interaction_class=interaction_class,
                    title=f"{title} (TCR alpha)",
                    xlabel=xlabel,
                    tcr_chain="alpha",
                    suffix="alpha",
                ),
                "beta": self._save_variant(
                    contacts=contacts,
                    output_root=output_root,
                    interaction_class=interaction_class,
                    title=f"{title} (TCR beta)",
                    xlabel=xlabel,
                    tcr_chain="beta",
                    suffix="beta",
                ),
            }

        return outputs
