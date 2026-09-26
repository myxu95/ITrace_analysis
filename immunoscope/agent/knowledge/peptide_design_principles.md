# Peptide Engineering Design Principles

## Overview

Peptide engineering for pMHC-TCR systems involves modifying the peptide sequence to:
1. **Enhance TCR recognition** (agonist design for vaccines, immunotherapy)
2. **Reduce TCR recognition** (antagonist design for autoimmunity)
3. **Improve MHC presentation** (vaccine design, epitope optimization)
4. **Modulate specificity** (cross-reactivity or specificity tuning)

Unlike TCR engineering, peptide engineering must balance **TCR recognition** with **MHC binding**, as the peptide must first be presented by MHC before TCR can recognize it.

## Peptide Position-Specific Rules

### Anchor Positions (P1, P2, PΩ-1, PΩ) - AVOID MUTATION

**Why avoid**:
- Buried in MHC groove
- Critical for peptide-MHC binding
- HLA allele-specific anchor motifs
- Mutation often abolishes MHC presentation

**Characteristics**:
- Low solvent exposure (<10%)
- High conservation across presented peptides
- Allele-specific preferences (e.g., HLA-A*02:01 prefers L, M at P2; V, L at PΩ)

**When mutation is acceptable**:
- Improving MHC binding (weak binders → strong binders)
- Conservative substitutions maintaining anchor properties
- Guided by MHC binding prediction tools

**Example**:
```
Original: P1=Y (suboptimal for HLA-A*02:01)
Improved: P1=L (preferred anchor for HLA-A*02:01)
Rationale: Increase MHC binding affinity for better presentation
```

### TCR Contact Positions (P4, P5, P6) - PRIMARY TARGET

**Why target**:
- Solvent-exposed (>70%)
- Direct TCR contact (CDR3 loops)
- High variability across peptides
- Determine TCR specificity

**Design strategies**:

#### 1. Enhance TCR Recognition (Agonist Design)

**Goal**: Increase TCR binding affinity for vaccine or immunotherapy

**Recommended mutations**:
- **Aromatic enhancement**: Introduce or optimize aromatic residues (F, Y, W) for π-π stacking with TCR CDR3
  - Example: P5 L→F (if TCR has aromatic residue at contact position)
  - Success rate: Moderate (40-60%)
  
- **Charge complementarity**: Introduce charged residues complementary to TCR
  - Example: P6 A→E (if TCR has R or K at contact position)
  - Success rate: Moderate (40-60%)
  
- **Hydrophobic optimization**: Enhance hydrophobic contacts
  - Example: P4 S→L (if TCR has hydrophobic pocket)
  - Success rate: Moderate (40-60%)

**Cautions**:
- Do not over-optimize: Super-agonists may cause T cell exhaustion
- Maintain MHC binding: Check that P4-P6 mutations don't disrupt MHC contact
- Consider immunogenicity: Non-self sequences may be more immunogenic (good for vaccines)

#### 2. Reduce TCR Recognition (Antagonist Design)

**Goal**: Decrease TCR binding for autoimmunity treatment or altered peptide ligands (APL)

**Recommended mutations**:
- **Disrupt key interactions**: Mutate hotspot positions to break critical contacts
  - Example: P5 F→A (remove π-π stacking)
  - Success rate: High (60-80%)
  
- **Introduce steric clashes**: Use bulky residues to prevent TCR binding
  - Example: P4 A→W (if space is limited)
  - Success rate: Moderate (40-60%)
  
- **Charge reversal**: Introduce repulsive electrostatic interactions
  - Example: P6 E→K (if TCR has K or R at contact position)
  - Success rate: High (60-80%)

**Cautions**:
- Maintain MHC binding: Antagonists must still be presented
- Avoid complete loss of binding: Partial agonists may be more useful than null peptides

#### 3. Modulate Specificity

**Goal**: Alter TCR recognition profile (cross-reactivity or specificity)

**Strategies**:
- **Conservative substitutions**: Maintain overall interaction type but tune strength
  - Example: P5 F→Y (both aromatic, but Y can H-bond)
  - Success rate: Moderate (40-60%)
  
- **Positional shifts**: Move key interaction to adjacent position
  - Example: P4 F→A, P5 A→F (shift aromatic contact)
  - Success rate: Low to moderate (30-50%)

### Flanking Positions (P3, P7) - SECONDARY TARGET

**Characteristics**:
- Partially exposed (30-50%)
- Contact both MHC and TCR (CDR1/2)
- Moderate variability

**Design strategies**:
- **Fine-tune TCR recognition**: Modulate CDR1/2 contacts
- **Adjust peptide conformation**: Influence P4-P6 presentation
- **Optimize MHC binding**: Contribute to peptide-MHC stability

**Recommended mutations**:
- Conservative substitutions maintaining size and charge
- Avoid drastic changes (may affect both MHC and TCR binding)

**Success rate**: Moderate (40-60%)

## MHC Binding Considerations

### Critical Principle

**Peptide must bind MHC before TCR can recognize it.**

Any peptide mutation must maintain sufficient MHC binding affinity (typically KD < 500 nM for effective presentation).

### MHC Binding Prediction

Use computational tools to predict MHC binding:
- **NetMHCpan**: Pan-specific MHC-I binding prediction
- **MHCflurry**: Machine learning-based prediction
- **IEDB**: Immune Epitope Database tools

**Workflow**:
1. Design peptide mutation
2. Predict MHC binding affinity
3. If predicted affinity is weak (KD > 500 nM), redesign
4. Prioritize mutations maintaining strong MHC binding

### HLA Allele-Specific Anchor Motifs

Different HLA alleles have different anchor preferences:

**HLA-A*02:01** (most common):
- P2: L, M, V, I, A (hydrophobic)
- PΩ: V, L, I, A (hydrophobic)

**HLA-A*01:01**:
- P2: T, S, D, E (small or charged)
- PΩ: Y, F (aromatic)

**HLA-A*03:01**:
- P2: L, M, V, I (hydrophobic)
- PΩ: K, R (basic)

**HLA-B*07:02**:
- P2: P (proline)
- PΩ: L, M, A (hydrophobic)

**Design rule**: Maintain or improve anchor residues according to HLA allele specificity.

## Peptide Design Strategies

### Strategy 1: Vaccine Design (Enhance Immunogenicity)

**Goal**: Create peptide vaccines that elicit strong T cell responses

**Workflow**:
1. **Identify weak epitopes**: Peptides with suboptimal MHC binding or TCR recognition
2. **Optimize MHC binding**: Improve anchor residues (P1, P2, PΩ)
3. **Enhance TCR recognition**: Optimize P4-P6 for stronger TCR contact
4. **Validate immunogenicity**: Test T cell activation in vitro/in vivo

**Example**:
```
Original peptide: SIINFEKL (OVA peptide, HLA-A*02:01)
Weak point: Suboptimal P2 anchor (I instead of L/M)
Improved: SLINFEKL (P2 I→L)
Result: 5-10x increase in MHC binding, enhanced T cell response
```

**Success rate**: High (60-80%) for MHC binding improvement

### Strategy 2: TCR Agonist Design

**Goal**: Design peptides that strongly activate specific TCR (for adoptive T cell therapy)

**Workflow**:
1. **Identify TCR hotspot contacts**: Use RRCS analysis to find key peptide positions
2. **Optimize hotspot interactions**: Enhance π-π stacking, salt bridges, hydrophobic contacts
3. **Maintain MHC binding**: Ensure anchor positions are not disrupted
4. **Test TCR activation**: Measure T cell activation (cytokine release, proliferation)

**Example**:
```
Original: P5=L (weak hydrophobic contact with TCR Y98)
Improved: P5=F (π-π stacking with TCR Y98)
Result: 3-5x increase in TCR binding, enhanced T cell activation
```

**Success rate**: Moderate (40-60%)

### Strategy 3: Altered Peptide Ligands (APL) for Autoimmunity

**Goal**: Design peptides that block pathogenic TCR without activating it (antagonists)

**Workflow**:
1. **Identify pathogenic TCR**: TCR causing autoimmune disease
2. **Design antagonist peptide**: Disrupt key TCR contacts at P4-P6
3. **Maintain MHC binding**: Ensure peptide is still presented
4. **Test antagonism**: Verify reduced T cell activation

**Example**:
```
Original: P5=F (strong π-π stacking with pathogenic TCR)
Antagonist: P5=A (remove π-π stacking)
Result: Peptide binds MHC but does not activate pathogenic TCR
```

**Success rate**: Moderate to high (50-70%)

### Strategy 4: Epitope Spreading (Multi-Epitope Vaccines)

**Goal**: Design multiple peptide variants to broaden T cell response

**Workflow**:
1. **Identify conserved epitope**: Peptide sequence conserved across pathogen strains
2. **Design variants**: Create peptides with different P4-P6 sequences
3. **Optimize each variant**: Enhance MHC binding and TCR recognition
4. **Test cross-reactivity**: Ensure variants activate different TCR clones

**Success rate**: Moderate (40-60%)

## Interaction-Specific Strategies for Peptide

### π-π Stacking

**Peptide side**:
- Introduce aromatic residues (F, Y, W) at P4-P6 if TCR has aromatic residues at contact positions
- F is most commonly used (smaller than W, no H-bond like Y)

**Example**:
```
TCR CDR3α Y98 contacts peptide P5
Current: P5=L (weak hydrophobic)
Improved: P5=F (π-π stacking)
Expected: 2-4x TCR binding improvement
```

### Salt Bridges

**Peptide side**:
- Introduce charged residues (E, D, K, R) at P4-P6 complementary to TCR
- E/D if TCR has K/R
- K/R if TCR has E/D

**Example**:
```
TCR CDR3β R155 contacts peptide P6
Current: P6=A (no interaction)
Improved: P6=E (salt bridge)
Expected: 2-3x TCR binding improvement
```

**Caution**: Charged residues at P4-P6 may also contact MHC helices - check for potential repulsion

### Hydrophobic Contacts

**Peptide side**:
- Optimize hydrophobic residues (L, I, V, M, F) at P4-P6 for TCR hydrophobic pockets
- Match size and shape complementarity

**Example**:
```
TCR has hydrophobic pocket at CDR3α
Current: P4=S (polar)
Improved: P4=L (hydrophobic)
Expected: 1.5-3x TCR binding improvement
```

## Common Pitfalls and How to Avoid Them

### Pitfall 1: Disrupting MHC Binding

**Problem**: Optimizing TCR contact but losing MHC presentation

**Solution**: 
- Always check MHC binding prediction after peptide mutation
- Prioritize mutations at P4-P6 (less likely to affect MHC)
- Avoid mutating anchor positions (P1, P2, PΩ) unless improving MHC binding

### Pitfall 2: Over-Optimization (Super-Agonists)

**Problem**: Too-strong TCR activation may cause T cell exhaustion or cytokine storm

**Solution**:
- Target moderate affinity improvement (3-10x, not 100x)
- Test T cell activation in vitro before in vivo
- Consider partial agonists for some applications

### Pitfall 3: Ignoring Immunogenicity

**Problem**: Non-self peptides may be immunogenic (good for vaccines, bad for therapeutics)

**Solution**:
- For vaccines: Non-self sequences are acceptable and may enhance response
- For therapeutics: Minimize deviation from self-sequences
- Use immunogenicity prediction tools

### Pitfall 4: Neglecting Proteasomal Processing

**Problem**: Peptide may not be generated by proteasome or transported by TAP

**Solution**:
- Check proteasomal cleavage sites (avoid introducing cleavage sites within epitope)
- Check TAP binding prediction (peptides must be transported to ER)
- Consider using minimal epitopes (8-11 aa) that bypass some processing steps

### Pitfall 5: Single-Position Focus

**Problem**: Optimizing one position while ignoring others

**Solution**:
- Consider multi-position optimization (but test single mutations first)
- Check for cooperative effects between positions
- Use MD simulation to assess overall peptide conformation

## Decision Framework for Peptide Design

### Step 1: Define Objective

- [ ] Enhance TCR recognition (agonist, vaccine)
- [ ] Reduce TCR recognition (antagonist, APL)
- [ ] Improve MHC presentation (vaccine, epitope optimization)
- [ ] Modulate specificity (cross-reactivity, specificity tuning)

### Step 2: Analyze MD Data

- [ ] Identify peptide hotspot positions (RRCS > 2.0 for peptide residues)
- [ ] Determine interaction types (H-bond, π-π, salt bridge, hydrophobic)
- [ ] Assess peptide flexibility (RMSF)
- [ ] Check peptide-MHC contacts (ensure stability)

### Step 3: Select Target Positions

- [ ] Prioritize P4-P6 (TCR contact, high tolerance)
- [ ] Consider P3, P7 (secondary targets)
- [ ] Avoid P1, P2, PΩ (MHC anchors) unless improving MHC binding

### Step 4: Design Mutations

- [ ] Use position-specific strategies (see above)
- [ ] Maintain or improve MHC binding (check prediction)
- [ ] Consider interaction type (π-π, salt bridge, hydrophobic)
- [ ] Design single mutations first

### Step 5: Predict Effects

- [ ] MHC binding prediction (NetMHCpan, MHCflurry)
- [ ] TCR binding prediction (if model available)
- [ ] Immunogenicity prediction (IEDB tools)

### Step 6: Prioritize and Test

- [ ] Rank mutations by predicted effect and confidence
- [ ] Test top 3-5 single mutations experimentally
- [ ] Measure MHC binding (competitive binding assay)
- [ ] Measure TCR activation (T cell assay)
- [ ] Consider combinations if single mutations successful

### Step 7: Iterate

- [ ] Analyze results (MHC binding, TCR activation, immunogenicity)
- [ ] Refine design based on experimental data
- [ ] Consider additional mutations or combinations

## Summary: Quick Reference Rules for Peptide Design

### DO:
✓ Focus on P4-P6 (TCR contact positions, high tolerance)
✓ Maintain MHC binding (check anchor positions P1, P2, PΩ)
✓ Use MHC binding prediction tools (NetMHCpan, MHCflurry)
✓ Consider HLA allele-specific anchor motifs
✓ Test single mutations before combinations
✓ Balance TCR recognition with MHC presentation
✓ Consider immunogenicity (good for vaccines, check for therapeutics)

### DON'T:
✗ Mutate anchor positions (P1, P2, PΩ) unless improving MHC binding
✗ Over-optimize TCR binding (super-agonists may cause exhaustion)
✗ Ignore MHC binding prediction (peptide must be presented)
✗ Neglect proteasomal processing and TAP transport
✗ Make multiple simultaneous mutations (unpredictable effects)
✗ Assume TCR optimization alone is sufficient (MHC binding is critical)

## Comparison: TCR vs Peptide Engineering

| Aspect | TCR Engineering | Peptide Engineering |
|--------|----------------|---------------------|
| **Primary target** | CDR3 loops | P4-P6 positions |
| **Avoid** | Framework | P1, P2, PΩ (anchors) |
| **Main constraint** | Structural stability | MHC binding |
| **Success rate** | 40-70% | 40-70% |
| **Validation** | TCR binding, T cell activation | MHC binding, TCR activation |
| **Application** | TCR therapy, TCR engineering | Vaccines, epitope optimization, APL |
| **Complexity** | High (large protein) | Moderate (short peptide) |
| **Predictability** | Moderate (structure-based) | Higher (MHC binding predictable) |
