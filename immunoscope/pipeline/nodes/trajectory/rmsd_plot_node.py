"""
RMSD Plot Node - Generate RMSD figure for preprocess quality reports.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class RMSDPlotNode(PipelineNode):
    """
    Generate RMSD curve plots for quality-report references.
    """

    def __init__(self,
                 output_filename: str = "rmsd.png",
                 output_subdir: str = "plots/rmsd",
                 title_prefix: str = "RMSD Assessment",
                 name: Optional[str] = None):
        super().__init__(name=name)
        self.output_filename = output_filename
        self.output_subdir = output_subdir
        self.title_prefix = title_prefix

    def validate_inputs(self, context: PipelineContext):
        rmsd_result = context.get_result('rmsd')
        if not rmsd_result:
            raise PipelineError(
                node_name=self.name,
                reason="rmsd result not found in context",
                context_state={'system_id': context.system_id},
            )

        output_file = rmsd_result.get('output_file')
        if not output_file or not Path(output_file).exists():
            raise PipelineError(
                node_name=self.name,
                reason="rmsd output file not found",
                context_state={
                    'system_id': context.system_id,
                    'output_file': output_file,
                },
            )

    def execute(self, context: PipelineContext) -> PipelineContext:
        self.logger.info(f"Executing {self.name} for system {context.system_id}")

        try:
            self.validate_inputs(context)
        except PipelineError as exc:
            context.add_error(str(exc))
            context.should_stop = True
            return context

        rmsd_result = context.get_result('rmsd')
        rmsd_file = Path(rmsd_result['output_file'])
        output_file = Path(context.get_subdir_path(self.output_subdir, self.output_filename))

        try:
            df = self._load_rmsd_dataframe(rmsd_file)
            time_ns = df.iloc[:, 0] / 1000.0
            rmsd_nm = df.iloc[:, 1]

            fig, ax = plt.subplots(figsize=(10, 5.5))
            ax.plot(time_ns, rmsd_nm, color="#1f4e79", linewidth=1.6, alpha=0.9)
            ax.axhline(rmsd_nm.mean(), color="#d97b29", linestyle="--", linewidth=1.2, alpha=0.9)

            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("RMSD (nm)")
            ax.set_title(f"{self.title_prefix} - {context.system_id}")
            ax.grid(True, alpha=0.25, linewidth=0.6)
            ax.set_ylim(bottom=0)

            stats_text = (
                f"Mean = {rmsd_nm.mean():.3f} nm\n"
                f"Std = {rmsd_nm.std():.3f} nm\n"
                f"Frames = {len(rmsd_nm)}"
            )
            ax.text(
                0.02,
                0.98,
                stats_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#cccccc"),
            )

            fig.tight_layout()
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close(fig)

            context.add_result(
                'rmsd_plot',
                {
                    'image_file': str(output_file),
                }
            )
            self.logger.info(f"RMSD plot generated: {output_file}")
        except Exception as exc:
            error_msg = f"RMSD plot generation failed: {exc}"
            context.add_error(error_msg)
            context.should_stop = True
            self.logger.exception(error_msg)

        return context

    def _load_rmsd_dataframe(self, rmsd_file: Path) -> pd.DataFrame:
        if rmsd_file.suffix.lower() == ".csv":
            return pd.read_csv(rmsd_file)

        rows = []
        with rmsd_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                rows.append((float(parts[0]), float(parts[1])))

        if not rows:
            raise ValueError(f"RMSD file contains no data: {rmsd_file}")

        return pd.DataFrame(rows, columns=["Time (ps)", "RMSD (nm)"])
