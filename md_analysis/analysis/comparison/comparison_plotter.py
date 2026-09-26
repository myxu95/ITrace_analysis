"""Plotting functions for comparison reports."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def write_quality_interface_comparison_plot(case_a_label: str, case_b_label: str, table: pd.DataFrame, output_file: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    metrics = [
        ("tail90_rmsd_variation_nm", "Tail-90% RMSD variation", "nm"),
        ("mean_bsa_angstrom2", "Mean BSA", "A^2"),
        ("mean_interface_ratio", "Mean interface ratio", ""),
    ]
    for axis, (metric, title, unit) in zip(axes, metrics):
        subset = table.loc[table["metric"] == metric]
        if subset.empty:
            _draw_empty(axis, title)
            continue
        row = subset.iloc[0]
        values = _valid_bar_values([row["case_a"], row["case_b"]])
        if not values:
            _draw_empty(axis, title)
            continue
        labels = [case_a_label, case_b_label]
        colors = ["#5d7d74", "#c08f5a"]
        axis.bar([labels[index] for index, _ in values], [value for _, value in values], color=[colors[index] for index, _ in values])
        axis.set_title(title, fontsize=11, weight="bold")
        if unit:
            axis.set_ylabel(unit)
        axis.grid(alpha=0.2, linestyle="--", axis="y")
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_flexibility_comparison_plot(case_a_label: str, case_b_label: str, frame: pd.DataFrame, output_file: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    frame = _filter_numeric_rows(frame, case_a_label, case_b_label)
    if frame.empty:
        _draw_empty(ax, "Region flexibility comparison")
    else:
        x = range(len(frame))
        width = 0.38
        ax.bar([i - width / 2 for i in x], frame[case_a_label].fillna(0.0), width=width, label=case_a_label, color="#5d7d74")
        ax.bar([i + width / 2 for i in x], frame[case_b_label].fillna(0.0), width=width, label=case_b_label, color="#c08f5a")
        ax.set_xticks(list(x))
        ax.set_xticklabels(frame["region_group"], rotation=30, ha="right")
        ax.set_ylabel("Mean RMSF (Å)")
        ax.set_title("Region-level flexibility comparison", fontsize=12, weight="bold")
        ax.grid(alpha=0.2, linestyle="--", axis="y")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_interaction_family_comparison_plot(case_a_label: str, case_b_label: str, frame: pd.DataFrame, output_file: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    frame = _filter_numeric_rows(frame, case_a_label, case_b_label)
    if frame.empty:
        _draw_empty(ax, "Interaction family comparison")
    else:
        x = range(len(frame))
        width = 0.38
        ax.bar([i - width / 2 for i in x], frame[case_a_label].fillna(0.0), width=width, label=case_a_label, color="#5d7d74")
        ax.bar([i + width / 2 for i in x], frame[case_b_label].fillna(0.0), width=width, label=case_b_label, color="#c08f5a")
        ax.set_xticks(list(x))
        ax.set_xticklabels(frame["family"], rotation=25, ha="right")
        ax.set_ylabel("Residue pairs")
        ax.set_title("Interaction family comparison", fontsize=12, weight="bold")
        ax.grid(alpha=0.2, linestyle="--", axis="y")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def write_rrcs_comparison_plot(case_a_label: str, case_b_label: str, frame: pd.DataFrame, output_file: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    frame = _filter_numeric_rows(frame, case_a_label, case_b_label)
    if frame.empty:
        _draw_empty(ax, "RRCS region comparison")
    else:
        x = range(len(frame))
        width = 0.38
        ax.bar([i - width / 2 for i in x], frame[case_a_label].fillna(0.0), width=width, label=case_a_label, color="#5d7d74")
        ax.bar([i + width / 2 for i in x], frame[case_b_label].fillna(0.0), width=width, label=case_b_label, color="#c08f5a")
        ax.set_xticks(list(x))
        ax.set_xticklabels(frame["interaction_class"], rotation=20, ha="right")
        ax.set_ylabel("Region RRCS sum")
        ax.set_title("RRCS region comparison", fontsize=12, weight="bold")
        ax.grid(alpha=0.2, linestyle="--", axis="y")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=220, facecolor="white")
    plt.close(fig)


def _draw_empty(axis, title: str) -> None:
    axis.set_title(title, fontsize=11, weight="bold")
    axis.text(0.5, 0.5, "No comparable data", ha="center", va="center", color="#666666", transform=axis.transAxes)
    axis.set_xticks([])
    axis.set_yticks([])


def _valid_bar_values(values: list[object]) -> list[tuple[int, float]]:
    valid: list[tuple[int, float]] = []
    for index, value in enumerate(values):
        numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if pd.notna(numeric):
            valid.append((index, float(numeric)))
    return valid


def _filter_numeric_rows(frame: pd.DataFrame, case_a_label: str, case_b_label: str) -> pd.DataFrame:
    if frame.empty or case_a_label not in frame.columns or case_b_label not in frame.columns:
        return pd.DataFrame()
    work = frame.copy()
    work[case_a_label] = pd.to_numeric(work[case_a_label], errors="coerce")
    work[case_b_label] = pd.to_numeric(work[case_b_label], errors="coerce")
    work = work.replace([np.inf, -np.inf], np.nan)
    return work[work[[case_a_label, case_b_label]].notna().any(axis=1)].reset_index(drop=True)
