#!/usr/bin/env python3
"""IMS Compare - two-system result comparison."""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

from immunoscope.core import PipelineContext
from immunoscope.analysis.comparison.cluster_representative_compare import build_cluster_representative_comparison
from immunoscope.pipeline.analysis_pipelines import SystemComparisonPipeline

# Add project root to path for scripts import
_project_root = Path(__file__).parent.parent.parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

from scripts.build_comparison_report_html import build_comparison_html_report
from scripts.build_enhanced_comparison_report import build_enhanced_comparison_html


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ims compare",
        description="Compare two analyzed ImmunoScope conditions and generate a comparison report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ims compare systems \
    --case-a output/md_standard \
    --case-b output/md_enhanced \
    --label-a "Standard sampling" \
    --label-b "Enhanced sampling" \
    --comparison-mode sampling \
    --comparison-scope same-system \
    --alignment-selection phla_core_ca \
    --comparison-context "Pre/post enhanced sampling" \
    -o output/compare_sampling
        """,
    )
    subparsers = parser.add_subparsers(dest="action", help="Comparison action")

    systems = subparsers.add_parser("systems", help="Compare two analyzed condition roots")
    systems.add_argument("--case-a", type=Path, required=True, help="Condition A result root or overview root")
    systems.add_argument("--case-b", type=Path, required=True, help="Condition B result root or overview root")
    systems.add_argument("--label-a", type=str, default="Condition A", help="Display label for condition A")
    systems.add_argument("--label-b", type=str, default="Condition B", help="Display label for condition B")
    systems.add_argument(
        "--comparison-mode",
        type=str,
        choices=["generic", "mutation", "sampling", "replicate"],
        default="generic",
        help="Comparison framing mode used for report wording",
    )
    systems.add_argument(
        "--comparison-context",
        type=str,
        default="",
        help="Optional free-text comparison context shown on the report cover",
    )
    systems.add_argument(
        "--comparison-scope",
        choices=["auto", "same-system", "cross-system"],
        default="auto",
        help="Comparison scope. Use same-system for Standard MD vs REST2 or replicate comparisons.",
    )
    systems.add_argument(
        "--alignment-selection",
        default="phla_core_ca",
        help="Declared fitting selection used or recommended before comparison, e.g. phla_core_ca or backbone.",
    )
    systems.add_argument(
        "--residue-mapping",
        choices=["auto", "identity", "region-only"],
        default="auto",
        help="Residue comparability policy for the report.",
    )
    systems.add_argument(
        "--alignment-report",
        type=Path,
        default=None,
        help="Optional alignment_report.json to embed into the comparison report.",
    )
    systems.add_argument("-o", "--output", type=Path, required=True, help="Comparison output directory")
    systems.add_argument(
        "--enhanced-report",
        action="store_true",
        help="Generate enhanced HTML report with hero metrics and priority-based layout",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def handle_systems(args, logger: logging.Logger) -> int:
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    context = PipelineContext(
        system_id=f"{args.label_a}_vs_{args.label_b}",
        topology="",
        trajectory_raw="",
        output_dir=str(output_dir),
    )
    pipeline = SystemComparisonPipeline(
        case_a_root=str(args.case_a),
        case_b_root=str(args.case_b),
        label_a=args.label_a,
        label_b=args.label_b,
        comparison_mode=args.comparison_mode,
        comparison_context=args.comparison_context,
        comparison_scope=args.comparison_scope,
        alignment_selection=args.alignment_selection,
        residue_mapping=args.residue_mapping,
    )

    logger.info("=" * 60)
    logger.info("IMS System Comparison")
    logger.info("=" * 60)
    logger.info("Condition A: %s", args.case_a)
    logger.info("Condition B: %s", args.case_b)
    logger.info("Comparison mode: %s", args.comparison_mode)
    logger.info("Comparison scope: %s", args.comparison_scope)
    logger.info("Alignment selection: %s", args.alignment_selection)
    if args.comparison_context:
        logger.info("Comparison context: %s", args.comparison_context)
    logger.info("Output: %s", output_dir)

    result_context = pipeline.execute(context)
    if result_context.has_errors():
        for error in result_context.errors:
            logger.error(error)
        return 1

    comparison_result = result_context.results["comparison"]
    cluster_representative_comparison = build_cluster_representative_comparison(
        case_a_root=args.case_a,
        case_b_root=args.case_b,
        case_a_label=args.label_a,
        case_b_label=args.label_b,
        output_dir=output_dir,
    )
    comparison_result["summary"]["cluster_representative_comparison"] = cluster_representative_comparison
    if cluster_representative_comparison.get("status") == "ready":
        comparison_result["artifacts"]["cluster_representative_comparison_json"] = cluster_representative_comparison.get("summary_json", "")
        comparison_result["artifacts"]["cluster_representative_rmsd_matrix_csv"] = cluster_representative_comparison.get("rmsd_matrix_csv", "")
        comparison_result["artifacts"]["cluster_representative_matches_csv"] = cluster_representative_comparison.get("matches_csv", "")
        comparison_result["artifacts"]["sampling_enriched_states_csv"] = cluster_representative_comparison.get("enriched_states_csv", "")
        comparison_result["artifacts"]["cluster_representative_rmsd_heatmap"] = cluster_representative_comparison.get("heatmap_png", "")
        comparison_result["artifacts"]["cluster_population_shift_plot"] = cluster_representative_comparison.get("population_shift_png", "")
    if args.alignment_report:
        alignment_report = json.loads(args.alignment_report.read_text(encoding="utf-8"))
        alignment_target = output_dir / "analysis" / "comparison" / "alignment_report.json"
        alignment_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.alignment_report, alignment_target)
        overlay_plot = alignment_report.get("structure_comparison", {}).get("overlay_plot")
        if overlay_plot and Path(overlay_plot).exists():
            overlay_target = output_dir / "analysis" / "comparison" / "initial_structure_overlay.png"
            shutil.copy2(overlay_plot, overlay_target)
            alignment_report["structure_comparison"]["overlay_plot"] = str(overlay_target.resolve())
            comparison_result["artifacts"]["initial_structure_overlay"] = str(overlay_target.resolve())
        for key, filename, artifact_key in (
            ("standard_pdb", "standard_frame0_protein.pdb", "initial_structure_standard_pdb"),
            ("rest2_aligned_pdb", "rest2_frame0_aligned_protein.pdb", "initial_structure_rest2_aligned_pdb"),
        ):
            pdb_path = alignment_report.get("structure_comparison", {}).get(key)
            if pdb_path and Path(pdb_path).exists():
                pdb_target = output_dir / "analysis" / "comparison" / filename
                shutil.copy2(pdb_path, pdb_target)
                alignment_report["structure_comparison"][key] = str(pdb_target.resolve())
                comparison_result["artifacts"][artifact_key] = str(pdb_target.resolve())
        comparison_result["summary"]["alignment_report"] = alignment_report
        comparison_result["artifacts"]["alignment_report_json"] = str(alignment_target.resolve())

    # Generate HTML report (enhanced or standard)
    if args.enhanced_report:
        html_path = build_enhanced_comparison_html(output_dir=output_dir, comparison_result=comparison_result)
        logger.info("Enhanced comparison report generated")
    else:
        html_path = build_comparison_html_report(output_dir=output_dir, comparison_result=comparison_result)
        logger.info("Comparison report generated")

    logger.info("HTML: %s", html_path)
    return 0


def main(argv=None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    if not args.action:
        parser.print_help()
        return 1

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    if args.action == "systems":
        return handle_systems(args, logger)

    logger.error("Unknown action: %s", args.action)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
