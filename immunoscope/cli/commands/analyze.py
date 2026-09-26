#!/usr/bin/env python3
"""IMS Analyze - grouped entry point for single-module analyses."""

from __future__ import annotations

from importlib import import_module
import sys


ANALYSIS_MODULES = {
    "preprocess": ("preprocess", "Process MD trajectories"),
    "quality": ("quality", "Assess trajectory quality"),
    "rmsd": ("rmsd", "Calculate RMSD"),
    "rmsf": ("rmsf", "Calculate RMSF"),
    "contact": ("contact", "Run coarse contact analysis"),
    "rrcs": ("rrcs", "Run residue-residue contact score analysis"),
    "identity": ("identity", "Annotate biological identity"),
    "bsa": ("bsa", "Calculate buried surface area"),
    "angle": ("angle", "Analyze docking angle"),
    "cluster": ("inter_cluster", "Run interface-aware state clustering"),
    "inter_cluster": ("inter_cluster", "Run interface-aware state clustering"),
    "landscape": ("landscape", "Build PCA/UMAP/TICA free energy landscapes"),
}


def print_usage() -> None:
    print(
        """
IMS Analyze - single-module analysis commands

Usage:
    ims analyze <module> [options]

Available modules:
    preprocess    Process MD trajectories
    quality       Assess trajectory quality
    rmsd          Calculate RMSD
    rmsf          Calculate RMSF
    contact       Run coarse contact analysis
    rrcs          Run residue-residue contact score analysis
    identity      Annotate biological identity
    bsa           Calculate buried surface area
    angle         Analyze docking angle
    cluster       Run interface-aware state clustering
    landscape     Build PCA/UMAP/TICA free energy landscapes

Examples:
    ims analyze rmsd -f processed.xtc -s md.tpr --selection backbone
    ims analyze rrcs --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/rrcs_demo
    ims analyze landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/landscape_demo --reducer tica

Legacy direct commands remain available, for example:
    ims rmsd ...
    ims landscape ...

For module help:
    ims analyze <module> --help
        """
    )


def _print_module_help(module_name: str, command_name: str) -> int:
    module = import_module(f"immunoscope.cli.commands.{module_name}")
    if not hasattr(module, "create_parser"):
        print(f"Module '{command_name}' does not expose parser help.")
        return 1
    parser = module.create_parser()
    parser.prog = f"ims analyze {command_name}"
    if parser.epilog:
        parser.epilog = parser.epilog.replace(f"ims {module_name}", f"ims analyze {command_name}")
        parser.epilog = parser.epilog.replace(f"ims {command_name}", f"ims analyze {command_name}")
    parser.print_help()
    return 0


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help", "help"}:
        print_usage()
        return 0

    command_name = args[0]
    if command_name not in ANALYSIS_MODULES:
        print(f"Error: Unknown analysis module '{command_name}'")
        print(f"\nAvailable modules: {', '.join(ANALYSIS_MODULES.keys())}")
        print("\nRun 'ims analyze --help' for more information.")
        return 1

    module_name, _description = ANALYSIS_MODULES[command_name]
    child_args = args[1:]
    if any(value in {"-h", "--help", "help"} for value in child_args):
        return _print_module_help(module_name, command_name)

    module = import_module(f"immunoscope.cli.commands.{module_name}")
    return module.main(child_args)


if __name__ == "__main__":
    raise SystemExit(main())
