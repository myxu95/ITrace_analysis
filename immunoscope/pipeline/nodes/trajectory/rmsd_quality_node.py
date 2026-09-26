"""
RMSD Quality Node - RMSD variation assessment and report generation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional

import pandas as pd

from ....core.base_node import PipelineNode
from ....core.context import PipelineContext
from ....core.exceptions import PipelineError


class RMSDQualityNode(PipelineNode):
    """
    Calculate variation from the RMSD time series and generate a report.
    """

    def __init__(self,
                 report_basename: str = "preprocess_quality_report",
                 name: Optional[str] = None):
        super().__init__(name=name)
        self.report_basename = report_basename

    def validate_inputs(self, context: PipelineContext):
        rmsd_result = context.get_result('rmsd')
        if not rmsd_result:
            raise PipelineError(
                node_name=self.name,
                reason="rmsd result not found in context",
                context_state={'system_id': context.system_id}
            )

        output_file = rmsd_result.get('output_file')
        if not output_file or not Path(output_file).exists():
            raise PipelineError(
                node_name=self.name,
                reason="rmsd output file not found",
                context_state={
                    'system_id': context.system_id,
                    'output_file': output_file,
                }
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

        try:
            times, rmsd_values = self._load_rmsd_series(rmsd_file)
            metrics = self._calculate_variation_metrics(times, rmsd_values)

            report_md = Path(context.get_quality_path(f"{self.report_basename}.md"))
            report_json = Path(context.get_quality_path(f"{self.report_basename}.json"))
            rmsd_plot = context.get_result('rmsd_plot')
            plot_file = None
            if rmsd_plot:
                plot_file = rmsd_plot.get('image_file')

            report_md.write_text(
                self._compose_markdown_report(
                    context=context,
                    metrics=metrics,
                    plot_file=plot_file,
                ),
                encoding='utf-8'
            )
            report_json.write_text(
                json.dumps(self._serialize_metrics(metrics, plot_file), indent=2, ensure_ascii=False),
                encoding='utf-8',
            )

            quality_result = {
                'metrics': self._serialize_metrics(metrics, plot_file),
                'reports': {
                    'markdown': str(report_md),
                    'json': str(report_json),
                },
            }
            context.add_result('preprocess_quality', quality_result)
            context.metadata['preprocess_quality'] = {
                'tail90_variation_nm': metrics['tail90_variation_nm'],
            }

            self.logger.info(
                f"RMSD variation assessment completed. "
                f"Tail90 variation: {metrics['tail90_variation_nm']:.3f} nm"
            )
        except Exception as exc:
            error_msg = f"RMSD variation assessment failed: {exc}"
            context.add_error(error_msg)
            context.should_stop = True
            self.logger.exception(error_msg)

        return context

    def _load_rmsd_series(self, rmsd_file: Path):
        if rmsd_file.suffix.lower() == '.csv':
            df = pd.read_csv(rmsd_file)
            return df.iloc[:, 0].to_numpy(), df.iloc[:, 1].to_numpy()

        times = []
        values = []
        with rmsd_file.open('r', encoding='utf-8') as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                times.append(float(parts[0]))
                values.append(float(parts[1]))

        if not times:
            raise ValueError(f"RMSD file contains no data: {rmsd_file}")

        return pd.Series(times).to_numpy(), pd.Series(values).to_numpy()

    def _calculate_variation_metrics(self, times, rmsd_values):
        n_frames = len(rmsd_values)
        if n_frames == 0:
            raise ValueError("RMSD series is empty")

        tail_start = min(n_frames - 1, max(0, math.floor(n_frames * 0.1)))
        tail_times = times[tail_start:]
        tail_values = rmsd_values[tail_start:]

        return {
            'n_frames': int(n_frames),
            'trim_fraction': 0.1,
            'trimmed_start_frame': int(tail_start),
            'trimmed_start_time_ps': float(times[tail_start]),
            'full_variation_nm': float(rmsd_values.max() - rmsd_values.min()),
            'tail90_variation_nm': float(tail_values.max() - tail_values.min()),
            'tail90_mean_rmsd_nm': float(tail_values.mean()),
            'tail90_std_rmsd_nm': float(tail_values.std()),
            'tail90_min_rmsd_nm': float(tail_values.min()),
            'tail90_max_rmsd_nm': float(tail_values.max()),
            'tail90_frame_count': int(len(tail_values)),
            'tail90_start_time_ps': float(tail_times[0]),
            'tail90_end_time_ps': float(tail_times[-1]),
        }

    def _serialize_metrics(self, metrics: dict, plot_file: Optional[str]) -> dict:
        serialized = {
            'plot_file': plot_file,
            'n_frames': int(metrics['n_frames']),
            'trim_fraction': float(metrics['trim_fraction']),
            'trimmed_start_frame': int(metrics['trimmed_start_frame']),
            'trimmed_start_time_ps': float(metrics['trimmed_start_time_ps']),
            'full_variation_nm': float(metrics['full_variation_nm']),
            'tail90_variation_nm': float(metrics['tail90_variation_nm']),
            'tail90_mean_rmsd_nm': float(metrics['tail90_mean_rmsd_nm']),
            'tail90_std_rmsd_nm': float(metrics['tail90_std_rmsd_nm']),
            'tail90_min_rmsd_nm': float(metrics['tail90_min_rmsd_nm']),
            'tail90_max_rmsd_nm': float(metrics['tail90_max_rmsd_nm']),
            'tail90_frame_count': int(metrics['tail90_frame_count']),
            'tail90_start_time_ps': float(metrics['tail90_start_time_ps']),
            'tail90_end_time_ps': float(metrics['tail90_end_time_ps']),
        }
        return serialized

    def _compose_markdown_report(self,
                                 context: PipelineContext,
                                 metrics: dict,
                                 plot_file: Optional[str]) -> str:
        lines = [
            f"# Preprocess Quality Report - {context.system_id}",
            "",
            f"**RMSD Full Variation**: {metrics['full_variation_nm']:.3f} nm",
            f"**RMSD Tail90 Variation**: {metrics['tail90_variation_nm']:.3f} nm",
            "",
        ]

        if plot_file:
            plot_path = Path(plot_file)
            report_dir = Path(context.get_quality_path("")).resolve()
            relative_plot = plot_path.resolve().relative_to(report_dir.parent.resolve())
            lines.extend([
                "## RMSD Plot",
                "",
                f"![RMSD Plot]({relative_plot.as_posix()})",
                "",
            ])

        lines.extend([
            "## RMSD Variation Summary",
            "",
            f"- Total frames: {metrics['n_frames']}",
            f"- Tail segment: last 90% (start frame {metrics['trimmed_start_frame']})",
            f"- Tail start time: {metrics['trimmed_start_time_ps']:.1f} ps",
            f"- Tail end time: {metrics['tail90_end_time_ps']:.1f} ps",
            f"- Full variation: {metrics['full_variation_nm']:.3f} nm",
            f"- Tail90 variation: {metrics['tail90_variation_nm']:.3f} nm",
            f"- Tail90 mean RMSD: {metrics['tail90_mean_rmsd_nm']:.3f} nm",
            f"- Tail90 std RMSD: {metrics['tail90_std_rmsd_nm']:.3f} nm",
            f"- Tail90 min RMSD: {metrics['tail90_min_rmsd_nm']:.3f} nm",
            f"- Tail90 max RMSD: {metrics['tail90_max_rmsd_nm']:.3f} nm",
            "",
        ])
        return "\n".join(lines)
