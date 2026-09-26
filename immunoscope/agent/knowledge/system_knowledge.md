# pMHC-TCR System Knowledge

## Complex Structure and Components

### Overview
The pMHC-TCR complex is a trimolecular assembly that mediates T cell recognition of antigenic peptides. The interaction is characterized by relatively weak binding affinity (KD ~1-100 μM) to allow rapid scanning and dissociation.

### Components

#### 1. Peptide
- **Length**: Typically 8-11 amino acids for MHC class I
- **Presentation**: Bound in extended conformation within MHC groove
- **Key positions**:
  - P1, P2 (N-terminus): Anchor residues, buried in MHC groove
  - P4-P6 (center): Primary TCR contact region, solvent-exposed
  - PΩ-1, PΩ (C-terminus): Anchor residues, buried in MHC groove
- **Variability**: Central positions (P4-P6) show high sequence diversity and determine TCR specificity

#### 2. MHC (HLA)
- **Structure**: α1 and α2 helices form peptide-binding groove
- **Function**: 
  - Presents peptide to TCR
  - Provides structural scaffold for peptide
  - Contributes to TCR binding surface
- **Polymorphism**: HLA alleles differ in groove residues, affecting peptide repertoire
- **TCR contact**: Mainly through α1/α2 helix residues flanking the peptide

#### 3. TCR (T Cell Receptor)
- **Structure**: αβ heterodimer with 6 CDR loops
- **CDR loops**:
  - **CDR1α/β**: Germline-encoded, contact MHC α2/α1 helices
  - **CDR2α/β**: Germline-encoded, contact MHC helices and peptide termini
  - **CDR3α/β**: Somatically recombined, highly variable, primary peptide contact
- **Variable domains**: Vα and Vβ provide framework for CDR loops
- **Constant domains**: Cα and Cβ mediate signaling (not in typical MD simulations)

### Interface Characteristics

#### Geometry
- **Buried Surface Area (BSA)**: 800-1200 Å² typical range
- **Contact residues**: 15-25 residue pairs across interface
- **Docking angle**: TCR typically crosses peptide-MHC at ~45-70° angle
- **Binding mode**: CDR3 loops centered over peptide, CDR1/2 contact MHC flanks

#### Interaction Network
- **Hydrogen bonds**: 8-15 H-bonds typical, distributed across CDRs
- **Salt bridges**: 2-5 charged pairs, often at periphery
- **Hydrophobic contacts**: Core interactions, especially CDR3-peptide
- **π-π stacking**: Aromatic residues in CDR3 with peptide aromatic/charged residues
- **Water-mediated**: Bridging waters stabilize interface

## Recognition Mechanism

### Specificity Determinants

#### CDR3 Loops (Primary)
- **CDR3α**: 
  - Contacts peptide center (P4-P6)
  - High sequence diversity (10-15 residues)
  - Determines fine specificity
  - Often contains aromatic residues for peptide recognition
  
- **CDR3β**:
  - Contacts peptide center and C-terminal region
  - High sequence diversity (12-18 residues)
  - Major contributor to binding energy
  - Can contact both peptide and MHC

#### CDR1/2 Loops (Secondary)
- **CDR1α/β**: Contact MHC α2/α1 helices, provide MHC restriction
- **CDR2α/β**: Contact MHC and peptide termini, modulate specificity
- **Germline-encoded**: Less variable, provide MHC-binding framework

### Binding Thermodynamics
- **Affinity range**: KD = 1-100 μM (relatively weak)
- **Kinetics**: Fast on-rate (10⁴-10⁶ M⁻¹s⁻¹), fast off-rate (0.01-1 s⁻¹)
- **Biological rationale**: Weak binding allows rapid scanning and serial triggering
- **Enthalpy/Entropy**: Often enthalpy-driven with unfavorable entropy (conformational restriction)

### Conformational Dynamics
- **CDR3 flexibility**: High RMSF (2-5 Å), allows induced fit
- **Peptide flexibility**: Central positions (P4-P6) show conformational sampling
- **Binding-induced changes**: CDR3 loops often undergo conformational selection
- **Allosteric effects**: Binding can propagate to constant domains (signaling)

## Functional Regions

### TCR Regions

#### Framework Regions (FR)
- **Location**: β-sheet scaffold supporting CDR loops
- **Function**: Maintain structural integrity, position CDR loops
- **Mutation tolerance**: Low - mutations may destabilize fold
- **Conservation**: High sequence conservation across TCRs

#### CDR Loops
- **CDR1/2**: 
  - Germline-encoded, moderate diversity
  - MHC-focused contacts
  - Mutation tolerance: Moderate
  - Design target: Secondary (for MHC cross-reactivity)

- **CDR3**:
  - Somatically recombined, extreme diversity
  - Peptide-focused contacts
  - Mutation tolerance: High (naturally variable)
  - Design target: Primary (for specificity tuning)

### Interface Regions

#### Core Interface
- **Definition**: Residues with >50% burial upon binding
- **Composition**: Hydrophobic and aromatic residues
- **Function**: Provide binding energy, shape complementarity
- **Mutation impact**: High - directly affects affinity

#### Peripheral Interface
- **Definition**: Residues with 20-50% burial
- **Composition**: Mix of polar and charged residues
- **Function**: Electrostatic steering, specificity tuning
- **Mutation impact**: Moderate - affects kinetics and specificity

#### Rim Region
- **Definition**: Residues near interface but not buried
- **Composition**: Charged and polar residues
- **Function**: Electrostatic pre-orientation, solubility
- **Mutation impact**: Low to moderate - affects association rate

## Key Concepts for Mutation Design

### Hotspot Residues
- **Definition**: Residues contributing >2 kcal/mol to binding energy
- **Identification**: RRCS analysis, alanine scanning, computational mutagenesis
- **Characteristics**: 
  - Often aromatic (Y, F, W) or charged (R, K, E, D)
  - Located in CDR3 loops
  - Form multiple interactions (H-bonds + hydrophobic)
- **Design relevance**: Prime targets for affinity enhancement

### Anchor Residues
- **Definition**: Residues maintaining structural integrity
- **Location**: Framework regions, CDR loop bases
- **Characteristics**: Conserved, buried, form internal H-bonds
- **Design relevance**: Avoid mutation - high risk of destabilization

### Specificity Residues
- **Definition**: Residues discriminating between similar peptides
- **Location**: CDR3 loops contacting variable peptide positions
- **Characteristics**: Variable across TCR repertoire, direct peptide contact
- **Design relevance**: Target for specificity tuning, cross-reactivity control

### Flexibility Regions
- **Definition**: Regions with high RMSF (>2 Å)
- **Location**: CDR3 loops, peptide center
- **Function**: Conformational adaptation, induced fit
- **Design relevance**: 
  - High flexibility: Tolerates mutations, allows optimization
  - Low flexibility: Risky to mutate, may be structurally constrained

## Biological Context

### T Cell Activation
- **Signal threshold**: Requires multiple pMHC-TCR engagements
- **Serial triggering**: TCR scans multiple pMHC, weak binding enables rapid dissociation
- **Kinetic proofreading**: Dwell time determines signaling outcome
- **Affinity window**: Too weak = no signal, too strong = reduced scanning efficiency

### Cross-Reactivity
- **Polyspecificity**: Single TCR recognizes ~10⁶ different peptides
- **Structural basis**: CDR3 flexibility allows accommodation of peptide variants
- **Design trade-off**: Increasing affinity may reduce cross-reactivity
- **Biological importance**: Enables broad pathogen coverage with limited TCR repertoire

### MHC Restriction
- **Definition**: TCR recognizes peptide only in context of specific MHC allele
- **Structural basis**: CDR1/2 loops contact MHC, CDR3 contacts peptide
- **Design implication**: Mutations affecting MHC contact may alter restriction
- **Clinical relevance**: TCR therapeutics must match patient HLA type

## Analysis Metrics and Interpretation

### RRCS (Residue-Residue Contact Strength)
- **Definition**: Weighted contact frequency × distance-based score
- **Interpretation**:
  - RRCS > 3.0: Strong hotspot, critical for binding
  - RRCS 1.5-3.0: Moderate contributor
  - RRCS < 1.5: Weak or transient contact
- **Design use**: Identify mutation targets, predict affinity impact

### BSA (Buried Surface Area)
- **Definition**: Surface area buried upon complex formation
- **Interpretation**:
  - BSA > 1000 Å²: Strong binding interface
  - BSA 800-1000 Å²: Typical TCR-pMHC
  - BSA < 800 Å²: Weak or partial interface
- **Design use**: Assess overall binding strength, guide interface optimization

### RMSF (Root Mean Square Fluctuation)
- **Definition**: Atomic positional fluctuation over trajectory
- **Interpretation**:
  - RMSF > 3 Å: High flexibility, induced fit region
  - RMSF 1-3 Å: Moderate flexibility, typical for CDR loops
  - RMSF < 1 Å: Rigid, structurally constrained
- **Design use**: Identify flexible regions tolerant to mutation

### Occupancy
- **Definition**: Fraction of trajectory with interaction present
- **Interpretation**:
  - Occupancy > 80%: Persistent, stable interaction
  - Occupancy 50-80%: Frequent but dynamic
  - Occupancy < 50%: Transient, may not be functionally important
- **Design use**: Distinguish stable from transient contacts

### Clustering
- **Definition**: Conformational states sampled during simulation
- **Interpretation**:
  - Single dominant cluster: Stable binding mode
  - Multiple clusters: Conformational heterogeneity, induced fit
  - Cluster transitions: Dynamic binding, flexibility
- **Design use**: Identify representative structures for design, assess binding stability
