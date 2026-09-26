"""Energy-landscape visualization and artifact writing."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .feature_matrix import FeatureMatrixResult
from .landscape_analyzer import LandscapeResult


class LandscapeVisualizer:
    """Generate static plots, interactive plots, and standard output files from analysis results."""

    def __init__(self, result: LandscapeResult):
        self.result = result

    def plot_2d_landscape(
        self,
        output_path: str,
        cmap: str = "itrace_fel",
        vmax: float | None = None,
        figsize: tuple[int, int] = (10, 8),
    ) -> None:
        fig, ax = plt.subplots(figsize=figsize)
        self._style_axis(ax)
        mesh_x, mesh_y = np.meshgrid(self.result.pc1_centers, self.result.pc2_centers)
        free_energy = self._display_free_energy(vmax=None)
        resolved_vmax = self._resolve_energy_vmax(free_energy, requested_vmax=vmax)
        color_map = self._landscape_cmap(cmap)
        contour = ax.contourf(
            mesh_x,
            mesh_y,
            free_energy.T,
            levels=self._energy_levels(resolved_vmax),
            cmap=color_map,
            vmax=resolved_vmax,
            extend="max",
        )
        self._add_colorbar(
            fig,
            ax,
            contour,
            label="Free energy",
            unit="kJ mol$^{-1}$",
            ticks=self._energy_ticks(resolved_vmax),
        )
        ax.set_xlabel(self._axis_label(0))
        ax.set_ylabel(self._axis_label(1))
        ax.set_title(f"{self.result.reducer.upper()} Free Energy Landscape")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    def plot_trajectory_overlay(
        self,
        output_path: str,
        stride: int = 10,
        alpha: float = 0.35,
        cmap: str = "itrace_fel",
        vmax: float | None = None,
    ) -> None:
        fig, ax = plt.subplots(figsize=(10, 8))
        self._style_axis(ax)
        mesh_x, mesh_y = np.meshgrid(self.result.pc1_centers, self.result.pc2_centers)
        free_energy = self._display_free_energy(vmax=None)
        resolved_vmax = self._resolve_energy_vmax(free_energy, requested_vmax=vmax)
        color_map = self._landscape_cmap(cmap)
        contour = ax.contourf(
            mesh_x,
            mesh_y,
            free_energy.T,
            levels=self._energy_levels(resolved_vmax),
            cmap=color_map,
            vmax=resolved_vmax,
            extend="max",
        )
        self._add_colorbar(
            fig,
            ax,
            contour,
            label="Free energy",
            unit="kJ mol$^{-1}$",
            ticks=self._energy_ticks(resolved_vmax),
        )

        coords = self.result.pc_coordinates
        ax.plot(coords[:, 0], coords[:, 1], color="#f7a531", linewidth=1.2, alpha=0.6)
        step = max(int(stride), 1)
        scatter = ax.scatter(
            coords[::step, 0],
            coords[::step, 1],
            c=np.arange(coords[::step].shape[0]),
            cmap="cividis",
            s=20,
            alpha=max(float(alpha), 0.55),
            linewidths=0.35,
            edgecolors="white",
        )
        self._add_colorbar(
            fig,
            ax,
            scatter,
            label="Trajectory order",
            orientation="horizontal",
            pad=0.11,
            fraction=0.045,
            aspect=38,
        )
        ax.set_xlabel(self._axis_label(0))
        ax.set_ylabel(self._axis_label(1))
        ax.set_title(f"Trajectory Overlay on {self.result.reducer.upper()} Landscape")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    def plot_loadings_heatmap(
        self,
        output_path: str,
        top_n: int = 10,
        figsize: tuple[int, int] = (10, 6),
    ) -> None:
        top_features = self._get_top_loading_feature_indices(top_n)
        if not top_features:
            raise ValueError("No feature loadings are available for plotting")

        selected_loadings = self.result.loadings[np.asarray(top_features, dtype=int), : min(2, self.result.loadings.shape[1])]
        selected_names = [self.result.feature_names[index] for index in top_features]

        fig, ax = plt.subplots(figsize=figsize)
        image = ax.imshow(selected_loadings, cmap="coolwarm", aspect="auto")
        label = "Loading" if self.result.reducer == "pca" else "Feature-embedding correlation"
        fig.colorbar(image, ax=ax, label=label)
        ax.set_xticks(range(selected_loadings.shape[1]))
        labels = self.result.coordinate_labels or [f"PC{i + 1}" for i in range(selected_loadings.shape[1])]
        ax.set_xticklabels(labels[: selected_loadings.shape[1]])
        ax.set_yticks(range(len(selected_names)))
        ax.set_yticklabels(selected_names)
        ax.set_title("Top Feature Loadings" if self.result.reducer == "pca" else "Top Feature-Embedding Correlations")

        for row_index in range(selected_loadings.shape[0]):
            for col_index in range(selected_loadings.shape[1]):
                ax.text(
                    col_index,
                    row_index,
                    f"{selected_loadings[row_index, col_index]:.2f}",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=8,
                )

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    def generate_interactive_html(
        self,
        output_path: str,
        include_trajectory: bool = True,
        vmax: float = 20.0,
    ) -> None:
        try:
            import plotly.graph_objects as go
        except ImportError as exc:  # pragma: no cover - runtime dependency
            raise ImportError("LandscapeVisualizer requires plotly to generate interactive plots.") from exc

        figure = go.Figure()
        figure.add_trace(
            go.Contour(
                z=self.result.free_energy.T,
                x=self.result.pc1_centers,
                y=self.result.pc2_centers,
                colorscale="Viridis",
                contours=dict(coloring="heatmap"),
                colorbar=dict(title="ΔG (kJ/mol)"),
                zmax=vmax,
                name="Free energy",
                hovertemplate=f"{self._coord_name(0)}=%{{x:.3f}}<br>{self._coord_name(1)}=%{{y:.3f}}<br>ΔG=%{{z:.3f}}<extra></extra>",
            )
        )
        if include_trajectory:
            figure.add_trace(
                go.Scatter(
                    x=self.result.pc_coordinates[:, 0],
                    y=self.result.pc_coordinates[:, 1],
                    mode="lines+markers",
                    marker=dict(size=4, color=np.arange(self.result.pc_coordinates.shape[0]), colorscale="Plasma"),
                    line=dict(color="rgba(247,165,49,0.7)", width=1.2),
                    name="Trajectory",
                    hovertemplate=f"{self._coord_name(0)}=%{{x:.3f}}<br>{self._coord_name(1)}=%{{y:.3f}}<extra></extra>",
                )
            )

        figure.update_layout(
            title=f"{self.result.reducer.upper()} Free Energy Landscape",
            xaxis_title=self._axis_label(0),
            yaxis_title=self._axis_label(1),
            template="plotly_white",
        )
        figure.write_html(output_path, include_plotlyjs=True)

    def write_report_assets(
        self,
        output_dir: str | Path,
        feature_result: FeatureMatrixResult | None = None,
        top_n: int = 10,
        include_interactive: bool = True,
    ) -> dict[str, str]:
        output_root = Path(output_dir)
        output_root.mkdir(parents=True, exist_ok=True)

        artifact_paths = {
            "feature_matrix": output_root / "feature_matrix.csv",
            "pca_coordinates": output_root / "pca_coordinates.csv",
            "loadings": output_root / "loadings.csv",
            "summary": output_root / "landscape_summary.json",
            "landscape_2d": output_root / "landscape_2d.png",
            "trajectory_overlay": output_root / "trajectory_overlay.png",
            "loadings_heatmap": output_root / "loadings_heatmap.png",
            "landscape_interactive": output_root / "landscape_interactive.html",
        }

        if feature_result is not None:
            feature_result.to_dataframe().to_csv(artifact_paths["feature_matrix"], index=False)
        self.result.to_pca_dataframe().to_csv(artifact_paths["pca_coordinates"], index=False)
        self.result.to_loadings_dataframe().to_csv(artifact_paths["loadings"], index=False)
        artifact_paths["summary"].write_text(
            json.dumps(
                {
                    "n_frames": int(self.result.pc_coordinates.shape[0]),
                    "n_features": int(len(self.result.feature_names)),
                    "n_components": int(self.result.n_components),
                    "reducer": self.result.reducer,
                    "coordinate_labels": self.result.coordinate_labels or [],
                    "explained_variance": [float(value) for value in self.result.explained_variance],
                    "reducer_metadata": self.result.reducer_metadata or {},
                    "top_contributors_pc1": [
                        {"feature": name, "loading": float(loading)}
                        for name, loading in self.result.get_top_contributors(0, top_n=top_n)
                    ],
                    "top_contributors_pc2": [
                        {"feature": name, "loading": float(loading)}
                        for name, loading in self.result.get_top_contributors(1, top_n=top_n)
                    ],
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        self.plot_2d_landscape(str(artifact_paths["landscape_2d"]))
        self.plot_trajectory_overlay(str(artifact_paths["trajectory_overlay"]))
        self.plot_loadings_heatmap(str(artifact_paths["loadings_heatmap"]), top_n=top_n)
        if include_interactive:
            self.generate_interactive_html(str(artifact_paths["landscape_interactive"]))

        output_mapping = {key: str(path) for key, path in artifact_paths.items()}
        if not include_interactive:
            output_mapping["landscape_interactive"] = ""
        return output_mapping

    def _get_top_loading_feature_indices(self, top_n: int) -> list[int]:
        n_components = min(2, self.result.loadings.shape[1])
        ranked_indices: list[int] = []
        for component_index in range(n_components):
            component = np.abs(self.result.loadings[:, component_index])
            ranked_indices.extend(np.argsort(component)[::-1][:top_n].tolist())

        ordered_unique_indices: list[int] = []
        for index in ranked_indices:
            if index not in ordered_unique_indices:
                ordered_unique_indices.append(index)
        return ordered_unique_indices[: max(top_n, 1)]

    def _coord_name(self, index: int) -> str:
        labels = self.result.coordinate_labels or []
        if index < len(labels):
            return labels[index]
        return f"PC{index + 1}"

    def _axis_label(self, index: int) -> str:
        label = self._coord_name(index)
        if self.result.reducer == "pca" and index < len(self.result.explained_variance):
            return f"{label} ({self.result.explained_variance[index] * 100:.1f}%)"
        return label

    @staticmethod
    def _style_axis(ax) -> None:
        ax.set_facecolor("#fbfcfd")
        ax.grid(True, color="#e6edf3", linewidth=0.7, alpha=0.75)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color("#cbd5e1")
            spine.set_linewidth(0.9)
        ax.tick_params(axis="both", colors="#334155", labelsize=9, length=3, width=0.8)
        ax.xaxis.label.set_color("#1f2937")
        ax.yaxis.label.set_color("#1f2937")
        ax.title.set_color("#111827")
        ax.title.set_fontweight("semibold")

    @staticmethod
    def _add_colorbar(
        fig,
        ax,
        mappable,
        label: str,
        unit: str = "",
        orientation: str = "vertical",
        ticks: list[float] | None = None,
        pad: float = 0.035,
        fraction: float = 0.045,
        aspect: int = 30,
    ):
        colorbar = fig.colorbar(
            mappable,
            ax=ax,
            orientation=orientation,
            fraction=fraction,
            pad=pad,
            aspect=aspect,
            ticks=ticks,
        )
        full_label = f"{label} ({unit})" if unit else label
        colorbar.set_label(full_label, fontsize=9, color="#334155", labelpad=8)
        colorbar.ax.tick_params(labelsize=8, colors="#475569", length=0)
        colorbar.outline.set_visible(False)
        if orientation == "vertical":
            colorbar.ax.yaxis.set_ticks_position("right")
            colorbar.ax.yaxis.set_label_position("right")
        else:
            colorbar.ax.xaxis.set_ticks_position("bottom")
            colorbar.ax.xaxis.set_label_position("bottom")
        colorbar.ax.set_facecolor("#fbfcfd")
        return colorbar

    @staticmethod
    def _energy_ticks(vmax: float | None) -> list[float] | None:
        if vmax is None or vmax <= 0:
            return None
        return [0.0, float(vmax) / 2.0, float(vmax)]

    @staticmethod
    def _energy_levels(vmax: float | None):
        if vmax is None or vmax <= 0:
            return 24
        return np.linspace(0.0, float(vmax), 25)

    @staticmethod
    def _resolve_energy_vmax(
        free_energy: np.ma.MaskedArray,
        requested_vmax: float | None = None,
    ) -> float:
        if requested_vmax is not None and requested_vmax > 0:
            return float(requested_vmax)

        values = np.ma.compressed(free_energy)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return 5.0

        percentile_value = float(np.percentile(values, 98))
        max_value = float(np.max(values))
        candidate = percentile_value if percentile_value > 0 else max_value
        candidate = max(candidate, 1.0)
        candidate = min(candidate, max_value if max_value > 0 else candidate)

        if candidate <= 2.0:
            return float(np.ceil(candidate * 4.0) / 4.0)
        if candidate <= 8.0:
            return float(np.ceil(candidate * 2.0) / 2.0)
        return float(np.ceil(candidate))

    def _display_free_energy(
        self,
        vmax: float | None = None,
        sigma: float = 1.0,
        support_sigma: float = 0.75,
        support_threshold: float = 0.015,
    ) -> np.ma.MaskedArray:
        """Build a report-friendly FEL surface without changing raw results.

        Histogram FEL is visually harsh for sparse trajectories. For plots we
        smooth the sampled probability field, recompute relative free energy,
        and mask regions with weak sampling support.
        """

        free_energy = np.asarray(self.result.free_energy, dtype=float)
        if self.result.histogram is None:
            return np.ma.array(free_energy, mask=~np.isfinite(free_energy))

        histogram = np.asarray(self.result.histogram, dtype=float)
        if histogram.shape != free_energy.shape or np.sum(histogram) <= 0:
            return np.ma.array(free_energy, mask=~np.isfinite(free_energy))

        try:
            from scipy.ndimage import gaussian_filter
        except ImportError:
            return self._masked_free_energy()

        probability = histogram / np.sum(histogram)
        smoothed_probability = gaussian_filter(probability, sigma=float(sigma), mode="nearest")
        smoothed_support = gaussian_filter((histogram > 0).astype(float), sigma=float(support_sigma), mode="nearest")
        low_support_mask = smoothed_support < float(support_threshold)

        positive = smoothed_probability[smoothed_probability > 0]
        if positive.size == 0:
            return self._masked_free_energy()

        probability_floor = float(np.percentile(positive, 5)) * 0.1
        safe_probability = np.where(smoothed_probability > 0, smoothed_probability, probability_floor)
        display_energy = -2.4942 * np.log(safe_probability)
        display_energy = display_energy - np.nanmin(display_energy)
        supported_values = display_energy[~low_support_mask & np.isfinite(display_energy)]
        if supported_values.size > 0:
            display_energy = np.where(
                low_support_mask,
                float(np.percentile(supported_values, 98)),
                display_energy,
            )
        if vmax is not None:
            display_energy = np.where(display_energy > vmax, float(vmax), display_energy)
        mask = ~np.isfinite(display_energy)
        return np.ma.array(display_energy, mask=mask)

    def _masked_free_energy(self) -> np.ma.MaskedArray:
        free_energy = np.asarray(self.result.free_energy, dtype=float)
        mask = ~np.isfinite(free_energy)
        if self.result.histogram is not None:
            histogram = np.asarray(self.result.histogram, dtype=float)
            if histogram.shape == free_energy.shape:
                mask = np.logical_or(mask, histogram <= 0)
        return np.ma.array(free_energy, mask=mask)

    @staticmethod
    def _landscape_cmap(name: str):
        if name == "itrace_fel":
            from matplotlib.colors import LinearSegmentedColormap

            color_map = LinearSegmentedColormap.from_list(
                "itrace_fel",
                [
                    "#1f5aa6",
                    "#2f7fc1",
                    "#6bb7d6",
                    "#d8edf2",
                    "#f7e0a3",
                    "#f2a15f",
                    "#c9443e",
                    "#7f1d1d",
                ],
                N=256,
            )
        else:
            color_map = plt.get_cmap(name)
        color_map = color_map.copy()
        color_map.set_bad(color="#7f1d1d", alpha=1.0)
        color_map.set_over("#6f1d1b")
        return color_map
