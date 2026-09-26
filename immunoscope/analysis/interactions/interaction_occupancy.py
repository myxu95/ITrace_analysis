"""
Generic interaction occupancy and stability analysis.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


class InteractionOccupancyAnalyzer:
    """Calculate pair-level stability and group-level contribution from annotated reports."""

    def analyze_report(self, report: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        if report.empty:
            return self._empty_pair_table(), self._empty_group_table()

        df = report.copy()
        for column, default in {
            "mhc_subregion": None,
            "frame_segments": None,
            "phla_region": None,
            "mhc_region": None,
            "tcr_chain": None,
            "tcr_region": None,
            "tcr_region_detailed": None,
            "interaction_class": None,
            "interaction_family": self._infer_interaction_family(report),
        }.items():
            if column not in df.columns:
                df[column] = default
        if "phla_resid" in df.columns:
            inferred_subregion = df.apply(self._infer_mhc_subregion_from_report, axis=1)
            missing_subregion = df["mhc_subregion"].isna() | (df["mhc_subregion"].astype(str).isin(["None", "nan"]))
            df.loc[missing_subregion, "mhc_subregion"] = inferred_subregion[missing_subregion]
        for col in [
            "contact_frequency",
            "contact_frames",
            "total_frames",
            "n_segments",
            "max_consecutive_frames",
            "max_consecutive_fraction",
        ]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        if "n_segments" not in df.columns:
            df["n_segments"] = 0
        if "max_consecutive_frames" not in df.columns:
            df["max_consecutive_frames"] = 0
        if "first_frame" not in df.columns:
            df["first_frame"] = None
        if "last_frame" not in df.columns:
            df["last_frame"] = None

        if "frame_segments" in df.columns:
            segment_summary = df["frame_segments"].apply(self._segment_summary_from_encoded)
            if "n_segments" in segment_summary.columns:
                missing = df["n_segments"].isna() | (pd.to_numeric(df["n_segments"], errors="coerce").fillna(0) == 0)
                df.loc[missing, "n_segments"] = segment_summary.loc[missing, "n_segments"]
            if "first_frame" in segment_summary.columns:
                missing = df["first_frame"].isna()
                df.loc[missing, "first_frame"] = segment_summary.loc[missing, "first_frame"]
            if "last_frame" in segment_summary.columns:
                missing = df["last_frame"].isna()
                df.loc[missing, "last_frame"] = segment_summary.loc[missing, "last_frame"]

        missing_max = df["max_consecutive_frames"].isna() | (pd.to_numeric(df["max_consecutive_frames"], errors="coerce").fillna(0) == 0)
        inferred_segments = pd.to_numeric(df["n_segments"], errors="coerce").fillna(0)
        contact_frames = pd.to_numeric(df["contact_frames"], errors="coerce").fillna(0)
        inferred_max = (contact_frames / inferred_segments.where(inferred_segments > 0, 1)).round().astype(int)
        df.loc[missing_max, "max_consecutive_frames"] = inferred_max[missing_max]

        if "max_consecutive_fraction" not in df.columns or df["max_consecutive_fraction"].isna().all():
            df["max_consecutive_fraction"] = df["max_consecutive_frames"] / df["total_frames"]
        df["max_consecutive_fraction"] = df["max_consecutive_fraction"].fillna(0.0)
        df["contact_frequency"] = df["contact_frequency"].fillna(0.0)
        df["n_segments"] = df["n_segments"].fillna(0).astype(int)
        df["max_consecutive_frames"] = df["max_consecutive_frames"].fillna(0).astype(int)
        df["contact_frames"] = df["contact_frames"].fillna(0).astype(int)

        df["stability_label"] = df.apply(self._classify_stability, axis=1)
        df["pair_label"] = df["phla_residue"].astype(str) + " ↔ " + df["tcr_residue"].astype(str)

        pair_table = df[
            [
                "interaction_family",
                "interaction_class",
                "pair_label",
                "phla_residue",
                "tcr_residue",
                "phla_region",
                "mhc_region",
                "mhc_subregion",
                "tcr_chain",
                "tcr_region",
                "tcr_region_detailed",
                "contact_frequency",
                "contact_frames",
                "total_frames",
                "n_segments",
                "max_consecutive_frames",
                "max_consecutive_fraction",
                "stability_label",
                "frame_segments",
            ]
        ].sort_values(
            by=["contact_frequency", "max_consecutive_fraction", "n_segments"],
            ascending=[False, False, True],
        ).reset_index(drop=True)

        group_summary = (
            df.groupby("interaction_class", as_index=False)
            .agg(
                n_pairs=("pair_label", "count"),
                mean_occupancy=("contact_frequency", "mean"),
                median_occupancy=("contact_frequency", "median"),
                mean_max_consecutive_fraction=("max_consecutive_fraction", "mean"),
                stable_pair_fraction=("stability_label", lambda s: float((s == "stable").mean()) if len(s) else 0.0),
                transient_pair_fraction=("stability_label", lambda s: float((s == "transient").mean()) if len(s) else 0.0),
            )
            .sort_values(by=["mean_occupancy", "n_pairs"], ascending=[False, False])
            .reset_index(drop=True)
        )

        return pair_table, group_summary

    @staticmethod
    def _infer_interaction_family(report: pd.DataFrame) -> str:
        if "interaction_family" in report.columns and not report["interaction_family"].dropna().empty:
            return str(report["interaction_family"].dropna().iloc[0])
        return "contact"

    def build_persistent_ranking(self, pair_table: pd.DataFrame) -> pd.DataFrame:
        if pair_table.empty:
            return pair_table.copy()
        return (
            pair_table.sort_values(
                by=["max_consecutive_fraction", "contact_frequency", "n_segments"],
                ascending=[False, False, True],
            )
            .reset_index(drop=True)
        )

    def build_residue_pair_occupancy_matrix(self, pair_table: pd.DataFrame) -> pd.DataFrame:
        if pair_table.empty:
            return pd.DataFrame()
        matrix = (
            pair_table.pivot_table(
                index="tcr_residue",
                columns="phla_residue",
                values="contact_frequency",
                aggfunc="max",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        matrix.index.name = "tcr_residue"
        return matrix

    def build_peptide_vs_mhc_summary(self, pair_table: pd.DataFrame) -> pd.DataFrame:
        if pair_table.empty:
            return self._empty_peptide_mhc_table()
        work = pair_table.copy()
        work["group"] = work["interaction_class"].map(
            {
                "peptide_tcr": "tcr_peptide",
                "hla_tcr": "tcr_mhc",
            }
        )
        work = work[work["group"].notna()].copy()
        if work.empty:
            return self._empty_peptide_mhc_table()
        return (
            work.groupby("group", as_index=False)
            .agg(
                n_pairs=("pair_label", "count"),
                mean_occupancy=("contact_frequency", "mean"),
                median_occupancy=("contact_frequency", "median"),
                mean_max_consecutive_fraction=("max_consecutive_fraction", "mean"),
                stable_pair_fraction=("stability_label", lambda s: float((s == "stable").mean()) if len(s) else 0.0),
                transient_pair_fraction=("stability_label", lambda s: float((s == "transient").mean()) if len(s) else 0.0),
            )
            .sort_values(by="mean_occupancy", ascending=False)
            .reset_index(drop=True)
        )

    def build_region_summary(self, pair_table: pd.DataFrame) -> pd.DataFrame:
        if pair_table.empty:
            return self._empty_region_table()
        work = pair_table.copy()
        work = work[work["tcr_region"].isin(["CDR1", "CDR2", "CDR3"])].copy()
        if work.empty:
            return self._empty_region_table()
        work["tcr_region_group"] = work["tcr_region"] + "_" + work["tcr_chain"].astype(str).str.lower()
        work["partner_region_group"] = work.apply(self._derive_partner_region_group, axis=1)
        work = work[work["partner_region_group"].isin(["peptide", "alpha1_helix", "alpha2_helix"])].copy()
        if work.empty:
            return self._empty_region_table()
        return (
            work.groupby(["tcr_region_group", "partner_region_group"], as_index=False)
            .agg(
                n_pairs=("pair_label", "count"),
                mean_occupancy=("contact_frequency", "mean"),
                median_occupancy=("contact_frequency", "median"),
                mean_max_consecutive_fraction=("max_consecutive_fraction", "mean"),
                stable_pair_fraction=("stability_label", lambda s: float((s == "stable").mean()) if len(s) else 0.0),
            )
            .sort_values(
                by=["tcr_region_group", "partner_region_group"],
                ascending=[True, True],
            )
            .reset_index(drop=True)
        )

    @staticmethod
    def write_summary_json(pair_table: pd.DataFrame, group_summary: pd.DataFrame, output_file: Path) -> None:
        summary = {
            "n_pairs": int(len(pair_table)),
            "n_stable_pairs": int((pair_table["stability_label"] == "stable").sum()) if not pair_table.empty else 0,
            "n_transient_pairs": int((pair_table["stability_label"] == "transient").sum()) if not pair_table.empty else 0,
            "interaction_classes": group_summary.to_dict(orient="records"),
        }
        output_file.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def build_all_outputs(self, report: pd.DataFrame) -> dict[str, pd.DataFrame]:
        pair_table, group_summary = self.analyze_report(report)
        persistent_ranking = self.build_persistent_ranking(pair_table)
        occupancy_matrix = self.build_residue_pair_occupancy_matrix(pair_table)
        peptide_mhc_summary = self.build_peptide_vs_mhc_summary(pair_table)
        region_summary = self.build_region_summary(pair_table)
        return {
            "pair_table": pair_table,
            "interaction_class_contribution": group_summary,
            "persistent_ranking": persistent_ranking,
            "occupancy_matrix": occupancy_matrix,
            "peptide_vs_mhc_summary": peptide_mhc_summary,
            "region_summary": region_summary,
        }

    @staticmethod
    def plot_outputs(pair_table: pd.DataFrame, group_summary: pd.DataFrame, output_dir: Path, family_name: str) -> dict[str, str]:
        output_dir.mkdir(parents=True, exist_ok=True)
        scatter_path = output_dir / f"{family_name}_occupancy_scatter.png"
        contribution_path = output_dir / f"{family_name}_interaction_class_contribution.png"

        fig, ax = plt.subplots(figsize=(8.4, 5.4))
        if pair_table.empty:
            ax.text(0.5, 0.5, "No occupancy records", ha="center", va="center", fontsize=14)
            ax.axis("off")
        else:
            palette = {
                "stable": "#2f5d50",
                "intermediate": "#c58a3d",
                "transient": "#c85a54",
            }
            for label, color in palette.items():
                sub = pair_table[pair_table["stability_label"] == label]
                if sub.empty:
                    continue
                ax.scatter(
                    sub["contact_frequency"],
                    sub["max_consecutive_fraction"],
                    s=42,
                    alpha=0.8,
                    color=color,
                    label=label,
                    edgecolors="white",
                    linewidths=0.6,
                )
            ax.set_xlabel("Occupancy / frequency")
            ax.set_ylabel("Longest consecutive fraction")
            ax.set_title(f"{family_name} occupancy stability")
            ax.set_xlim(0, 1.05)
            ax.set_ylim(0, 1.05)
            ax.grid(alpha=0.18, linestyle="--")
            ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(scatter_path, dpi=220, facecolor="white")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.4, 5.2))
        if group_summary.empty:
            ax.text(0.5, 0.5, "No group contribution data", ha="center", va="center", fontsize=14)
            ax.axis("off")
        else:
            plot_df = group_summary.sort_values(by="mean_occupancy", ascending=True)
            ax.barh(plot_df["interaction_class"], plot_df["mean_occupancy"], color="#5b8f7d")
            for idx, row in plot_df.iterrows():
                ax.text(
                    float(row["mean_occupancy"]) + 0.01,
                    idx,
                    f"stable {row['stable_pair_fraction']:.0%}",
                    va="center",
                    fontsize=10,
                    color="#2d3a35",
                )
            ax.set_xlim(0, 1.05)
            ax.set_xlabel("Mean occupancy")
            ax.set_title(f"{family_name} class contribution")
            ax.grid(axis="x", alpha=0.18, linestyle="--")
            ax.set_axisbelow(True)
        fig.tight_layout()
        fig.savefig(contribution_path, dpi=220, facecolor="white")
        plt.close(fig)

        return {
            "scatter_plot": str(scatter_path),
            "contribution_plot": str(contribution_path),
        }

    def plot_timeline(self, pair_table: pd.DataFrame, output_dir: Path, family_name: str, top_n: int = 8) -> str:
        output_dir.mkdir(parents=True, exist_ok=True)
        timeline_path = output_dir / f"{family_name}_persistent_timeline.png"
        fig, ax = plt.subplots(figsize=(10.2, 4.8))
        if pair_table.empty:
            ax.text(0.5, 0.5, "No interaction timeline records", ha="center", va="center", fontsize=14)
            ax.axis("off")
        else:
            top = self.build_persistent_ranking(pair_table).head(top_n).copy()
            plotted = False
            for idx, (_, row) in enumerate(top.iterrows()):
                segments = self._decode_segments(row.get("frame_segments"))
                if not segments:
                    continue
                bars = [(start, max(end - start, 1)) for start, end in segments]
                ax.broken_barh(bars, (idx - 0.35, 0.7), facecolors="#5b8f7d", edgecolors="#2f5d50", linewidth=0.8)
                plotted = True
            if plotted:
                ax.set_yticks(range(len(top)))
                ax.set_yticklabels(top["pair_label"].tolist(), fontsize=9)
                ax.set_xlabel("Frame index")
                ax.set_title(f"{family_name} persistent interaction timeline")
                ax.grid(axis="x", alpha=0.18, linestyle="--")
                ax.invert_yaxis()
            else:
                ax.text(0.5, 0.5, "No decoded timeline segments", ha="center", va="center", fontsize=14)
                ax.axis("off")
        fig.tight_layout()
        fig.savefig(timeline_path, dpi=220, facecolor="white")
        plt.close(fig)
        return str(timeline_path)

    @staticmethod
    def _classify_stability(row: pd.Series) -> str:
        occupancy = float(row.get("contact_frequency", 0.0) or 0.0)
        longest = float(row.get("max_consecutive_fraction", 0.0) or 0.0)
        segments = int(row.get("n_segments", 0) or 0)
        if occupancy >= 0.30 and longest >= 0.10:
            return "stable"
        if occupancy < 0.10 and longest < 0.03 and segments >= 2:
            return "transient"
        return "intermediate"

    @staticmethod
    def _derive_partner_region_group(row: pd.Series) -> str | None:
        if str(row.get("interaction_class")) == "peptide_tcr":
            return "peptide"
        subregion = row.get("mhc_subregion")
        if subregion in {"alpha1_helix", "alpha2_helix"}:
            return str(subregion)
        return None

    @staticmethod
    def _infer_mhc_subregion_from_report(row: pd.Series) -> str | None:
        if str(row.get("interaction_class")) == "peptide_tcr":
            return "peptide"
        if str(row.get("phla_region")) != "HLA_alpha":
            return None
        try:
            resid = int(row.get("phla_resid"))
        except (TypeError, ValueError):
            return None
        if 50 <= resid <= 86:
            return "alpha1_helix"
        if 140 <= resid <= 176:
            return "alpha2_helix"
        return "non_groove"

    @staticmethod
    def _decode_segments(raw_segments: Any) -> list[list[int]]:
        if raw_segments is None or (isinstance(raw_segments, float) and pd.isna(raw_segments)):
            return []
        if isinstance(raw_segments, list):
            return [[int(item[0]), int(item[1])] for item in raw_segments if isinstance(item, (list, tuple)) and len(item) == 2]
        try:
            parsed = json.loads(str(raw_segments))
        except (TypeError, ValueError, json.JSONDecodeError):
            return []
        if not isinstance(parsed, list):
            return []
        cleaned: list[list[int]] = []
        for item in parsed:
            if isinstance(item, list) and len(item) == 2:
                cleaned.append([int(item[0]), int(item[1])])
        return cleaned

    def _segment_summary_from_encoded(self, raw_segments: Any) -> pd.Series:
        segments = self._decode_segments(raw_segments)
        if not segments:
            return pd.Series(
                {
                    "n_segments": 0,
                    "first_frame": None,
                    "last_frame": None,
                }
            )
        return pd.Series(
            {
                "n_segments": len(segments),
                "first_frame": int(segments[0][0]),
                "last_frame": int(segments[-1][1]),
            }
        )

    @staticmethod
    def _empty_pair_table() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "interaction_family",
                "interaction_class",
                "pair_label",
                "phla_residue",
                "tcr_residue",
                "phla_region",
                "mhc_region",
                "mhc_subregion",
                "tcr_chain",
                "tcr_region",
                "tcr_region_detailed",
                "contact_frequency",
                "contact_frames",
                "total_frames",
                "n_segments",
                "max_consecutive_frames",
                "max_consecutive_fraction",
                "stability_label",
                "frame_segments",
            ]
        )

    @staticmethod
    def _empty_group_table() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "interaction_class",
                "n_pairs",
                "mean_occupancy",
                "median_occupancy",
                "mean_max_consecutive_fraction",
                "stable_pair_fraction",
                "transient_pair_fraction",
            ]
        )

    @staticmethod
    def _empty_peptide_mhc_table() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "group",
                "n_pairs",
                "mean_occupancy",
                "median_occupancy",
                "mean_max_consecutive_fraction",
                "stable_pair_fraction",
                "transient_pair_fraction",
            ]
        )

    @staticmethod
    def _empty_region_table() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "tcr_region_group",
                "partner_region_group",
                "n_pairs",
                "mean_occupancy",
                "median_occupancy",
                "mean_max_consecutive_fraction",
                "stable_pair_fraction",
            ]
        )
