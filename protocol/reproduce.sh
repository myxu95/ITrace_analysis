#!/usr/bin/env bash
# ============================================================================
# ITrace — reproduce one complex from crystal structure to analysis-ready MD.
#   Build  ->  Equilibrate  ->  3 x 200 ns production  ->  post-process
#
# Requirements: GROMACS 2026.x, the CHARMM36 (Jul 2021) force field installed,
# and the four .mdp files from this bundle (em/nvt/npt/md) in the working dir.
# Production is GPU-heavy (~200 ns x 3); run each replica on a GPU node.
# This is a template — adapt PDB id, box size and the cleanup step to your case.
# ============================================================================
set -euo pipefail

PDB=1ao7                     # crystal structure to reproduce

# --- 0. Clean the crystal: keep the biological assembly, drop crystallographic
#        waters / ligands, and add any missing heavy atoms (e.g. with PDBFixer).
#        Produces clean.pdb with chains A=MHC B=b2m C=peptide D=TCRa E=TCRb.
#        (See the tutorial; result assumed as clean.pdb below.)

# --- 1. Topology: CHARMM36 / TIP3P ---
gmx pdb2gmx -f clean.pdb -o proc.gro -p topol.top -water tip3p -ff charmm36-jul2021 -ignh

# --- 2. Box, solvate, neutralise (0.15 M NaCl) ---
gmx editconf -f proc.gro -o box.gro -c -d 1.0 -bt cubic
gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp   -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                      -neutral -conc 0.15 -pname NA -nname CL

# --- 3. Energy minimisation ---
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun  -deffnm em

# --- 4. Equilibration: NVT (100 ps) then NPT (1 ns) ---
gmx grompp -f nvt.mdp -c em.gro  -r em.gro  -p topol.top -o nvt.tpr
gmx mdrun  -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun  -deffnm npt

# --- 5. Production: 3 independent replicas (independent initial velocities) ---
for r in 1 2 3; do
  gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_run${r}.tpr
  gmx mdrun  -deffnm md_run${r}                       # add -gpu_id / -nb gpu on a GPU node
done

# --- 6. Post-process each replica into the ITrace served format ------------
#        PBC: whole -> nojump -> fit(rot+trans on backbone); then strip water/
#        ions to a protein-only trajectory + a topology PDB.
for r in 1 2 3; do
  gmx trjconv -s md_run${r}.tpr -f md_run${r}.xtc -o tmp_whole.xtc  -pbc whole  <<< "System"
  gmx trjconv -s md_run${r}.tpr -f tmp_whole.xtc  -o tmp_nojump.xtc -pbc nojump <<< "System"
  gmx trjconv -s md_run${r}.tpr -f tmp_nojump.xtc -o md_run${r}_processed.xtc \
              -fit rot+trans <<< $'Backbone\nProtein'
  gmx trjconv -s md_run${r}.tpr -f md_run${r}_processed.xtc -o run${r}_topology.pdb \
              -dump 0 <<< "Protein"
  rm -f tmp_whole.xtc tmp_nojump.xtc
done

echo "Done. md_run{1,2,3}_processed.xtc + run{1,2,3}_topology.pdb are analysis-ready."
