#!/usr/bin/env python3
"""
IMS - ImmunoScope Command Line Interface

Main entry point for all ImmunoScope CLI commands.

Author: ImmunoScope Development Team
Date: 2026-03-15
"""

import argparse
import sys
from importlib import import_module


RECOMMENDED_SUBCOMMANDS = {
    'run': 'Run guided analysis workflows',
    'analyze': 'Run a single analysis module',
    'batch': 'Run batch processing workflows',
    'report': 'Generate or serve analysis reports',
    'agent': 'Interactive agent for MD analysis',
    'agent-web': 'Run the web agent interface',
    'compare': 'Compare analyzed systems',
    # 2026-05-27: `recommend` 子命令已退役。一次性 batch 推荐这种形态
    # 已被对话式的 Design Copilot（agent-web 路径）替代——`web/routers/
    # design.py` 顶部那行注释早就说了 "Replaces the old batch
    # recommendation engine"，但当时 CLI 入口没真的撤掉。这次清掉。
    'web': 'Run the ImmunoScope web interface',
    'pdb': 'PDB file processing utilities',
}

LEGACY_DIRECT_SUBCOMMANDS = {
    'preprocess': 'Process MD trajectories (PBC correction)',
    'quality': 'Quality assessment and validation',
    'rmsd': 'RMSD calculation and analysis',
    'rmsf': 'RMSF calculation and analysis',
    'contact': 'Contact analysis',
    'rrcs': 'Residue-residue contact score analysis',
    'identity': 'Biological identity annotation',
    'bsa': 'Buried surface area analysis',
    'residue_sasa': 'Per-residue bound/unbound SASA analysis',
    'dssp': 'Per-residue secondary structure (DSSP)',
    'chi_dihedrals': 'Sidechain chi1/chi2 dihedral entropy',
    'inter_cluster': 'Interface-aware state clustering',
    'angle': 'Docking angle analysis',
    'landscape': 'Energy landscape analysis',
}

SUBCOMMANDS = {
    **RECOMMENDED_SUBCOMMANDS,
    **LEGACY_DIRECT_SUBCOMMANDS,
}


def print_usage():
    """Print usage information"""
    print("""
IMS - ImmunoScope MD Analysis Toolkit

Usage:
    ims <subcommand> [options]

Recommended subcommands:
    run           Run guided analysis workflows
    analyze       Run a single analysis module
    batch         Run batch processing workflows
    report        Generate or serve analysis reports
    agent         Interactive agent for MD analysis
    agent-web     Run the web agent interface
    compare       Compare analyzed systems
    web           Run the ImmunoScope web interface
    pdb           PDB file processing utilities

Legacy direct analysis commands:
    preprocess    Process MD trajectories (PBC correction)
    quality       Quality assessment and validation
    rmsd          RMSD calculation and analysis
    rmsf          RMSF calculation and analysis
    contact       Contact analysis
    rrcs          Residue-residue contact score analysis
    identity      Biological identity annotation
    bsa           Buried surface area analysis
    residue_sasa  Per-residue bound/unbound SASA analysis
    dssp          Per-residue secondary structure (DSSP)
    chi_dihedrals Sidechain chi1/chi2 dihedral entropy
    inter_cluster Interface-aware state clustering
    angle         Docking angle analysis
    landscape     Energy landscape analysis

Examples:
    ims run single --structure raw.pdb --topology md.tpr --trajectory md.xtc -o output/job_001
    ims preprocess -f md.xtc -s md.tpr -o processed.xtc
    ims quality -f processed.xtc -s md.tpr
    ims analyze rmsd -f processed.xtc -s md.tpr --selection backbone
    ims analyze rrcs --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/rrcs_demo
    ims analyze identity --structure md_processed_converted.pdb -o ./output/identity_demo
    ims analyze bsa -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./bsa
    ims analyze cluster -f md_processed.xtc -s md.tpr --structure md_processed_converted.pdb -o ./output/interface_cluster
    ims analyze landscape --structure md_processed_converted.pdb --topology md.tpr --trajectory md_processed.xtc -o ./output/landscape_demo
    ims compare systems --case-a output/case_A --case-b output/case_B -o output/compare_ab
    ims web --host 127.0.0.1 --port 8890
    ims agent-web --port 8891       # browser UI for the Agent, SSH-forward friendly
    ims report interaction --base-dir output/interaction_case_1OGA --system-id 1OGA_sd_run2
    ims report serve --overview-dir output/1OGA_demo_bundle/report/interaction_case_1OGA/overview
    ims agent                # interactive MD analysis agent (replaces `assistant`)
    ims pdb download 1ao7 2ckb
    ims batch preprocess /data/md/ --workers 4

For help on specific subcommand:
    ims analyze --help
    ims analyze landscape --help
    ims agent --help
    """)


def main():
    """Main entry point for IMS CLI"""

    # If no arguments, print usage
    if len(sys.argv) < 2:
        print_usage()
        return 0

    # Check if subcommand is valid
    subcommand = sys.argv[1]

    # Handle help
    if subcommand in ['-h', '--help', 'help']:
        print_usage()
        return 0

    # Handle version
    if subcommand in ['-v', '--version', 'version']:
        from immunoscope import __version__
        print(f"IMS (ImmunoScope) version {__version__}")
        return 0

    # Validate subcommand
    if subcommand not in SUBCOMMANDS:
        print(f"Error: Unknown subcommand '{subcommand}'")
        print(f"\nAvailable subcommands: {', '.join(SUBCOMMANDS.keys())}")
        print(f"\nRun 'ims --help' for more information.")
        return 1

    # Import and run subcommand
    try:
        module = import_module(f'immunoscope.cli.commands.{subcommand}')
        # Pass remaining arguments to subcommand
        return module.main(sys.argv[2:])
    except ImportError as e:
        print(f"Error loading subcommand '{subcommand}': {e}")
        import traceback
        traceback.print_exc()
        return 1
    except Exception as e:
        print(f"Error running subcommand '{subcommand}': {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
