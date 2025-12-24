# Simulation Pipeline Description

Our simulation procedure consists of four major phases: (i) simulating species trees, (ii) simulating gene trees with duplication/loss and incomplete lineage sorting, (iii) simulating DNA sequences and estimating gene trees, and (iv) simulating protein sequences and estimating gene trees. The selection of parameters in all phases is based on the simulation study by Molloy and Warnow (2020), which uses parameters derived from a fungal dataset presented by Rasmussen and Kellis (2011).

## Phase 1: Simulating Species Trees

Species trees are simulated using **SimPhy** version 1.0.2+ with a birth-death process using parameters from Molloy and Warnow (2020):

- **Tree height**: 1,800,000,337.5 years
- **Speciation rate**: λ = 1.8 × 10⁻⁹ events/year
- **Extinction rate**: μ = 0 events/year
- **Number of leaves**: 12 or 20 (configurable via command-line)

SimPhy produces **ultrametric species trees**, where all leaves are equidistant from the root.

### SimPhy Command for Species Trees

```bash
simphy_lnx64 \
  -rs <num_replicates> \
  -sb f:<speciation_rate> \
  -sd f:<extinction_rate> \
  -sl f:<num_leaves> \
  -stl f:<tree_height> \
  -sp f:100000 \
  -o <output_dir>
```

Trees are extracted from SimPhy output (`s_tree.trees`) and saved in simple Nexus (.nex) format for use in Phase 2.

## Phase 2: Simulating Gene Trees

Gene trees are simulated using **SimPhy** version 1.0.2+ to model:
- **Duplication and loss (DL) events** along species tree branches
- **Incomplete lineage sorting (ILS)** through coalescent simulation

### Key Parameters

- **Substitution rate**: 4 × 10⁻¹⁰ substitutions/site/year
- **Duplication/loss rates**: {10⁻¹⁰, 2×10⁻¹⁰, 5×10⁻¹⁰} events/year
- **Population sizes**: {10⁷, 5×10⁷}

SimPhy uses these rates directly (in years) and outputs gene tree branch lengths already scaled to **substitutions per site**.

### SimPhy Command

```bash
simphy_lnx64 \
  -sr <species_tree.nex> \
  -rg 1 \
  -rl F:1 \
  -si F:1 \
  -sp F:<population_size> \
  -su F:0.0000000004 \
  -lb F:<dl_rate> \
  -ld F:lb \
  -ll 3 \
  -hg LN:1.5,1 \
  -oc 1 \
  -o <output_folder> \
  -v 0 \
  -cs <unique_seed>
```

**Parameters:**
- `-sr`: Species tree in Nexus format
- `-rg 1`: Generate 1 gene tree per locus
- `-rl F:1`: Generate 1 locus
- `-si F:1`: Generate 1 individual per species
- `-sp F:<pop_size>`: Effective population size (10⁷ or 5×10⁷)
- `-su F:0.0000000004`: Substitution rate (4×10⁻¹⁰ substitutions/site/year)
- `-lb F:<dl_rate>`: Birth (duplication) rate (e.g., 1e-10, 2e-10, or 5e-10 events/year)
- `-ld F:lb`: Death (loss) rate = birth rate
- `-ll 3`: Minimum number of lineages
- `-hg LN:1.5,1`: Gene tree height follows lognormal distribution
- `-oc 1`: Output coalescent trees
- `-cs <unique_seed>`: Random seed (varies by configuration and replicate for independence)

### Parameter Combinations

Following Molloy and Warnow (2020), we use:

**Duplication/Loss Rates** (3 levels):
- Low: 10⁻¹⁰ events/year (biological data from *Saccharomycetales* fungi)
- Medium: 2×10⁻¹⁰ events/year
- High: 5×10⁻¹⁰ events/year

**Population Sizes** (2 levels):
- Low ILS: 10⁷ (corresponding to low incomplete lineage sorting)
- Medium ILS: 5×10⁷ (corresponding to medium incomplete lineage sorting)

This yields **6 parameter combinations** per leaf-set size, and **12 total configurations** for 2 leaf-set sizes.

### Critical Implementation Details

1. **Direct rate usage**: SimPhy uses rates directly in years (e.g., 1e-10, not scaled)
2. **Independent seeds**: Each configuration and replicate uses a unique seed to ensure independent stochastic events
3. **Branch lengths in substitutions**: SimPhy outputs gene tree branch lengths in **substitutions per site** (not time units)
4. **Output**: SimPhy produces three tree files per simulation:
   - `s_tree.trees`: Species tree
   - `l_trees.trees`: Locus tree (after duplication/loss)
   - `g_trees1.trees`: Gene tree (after ILS coalescence)

## Phase 3: Simulating DNA Sequences and Estimating Gene Trees

DNA sequences are simulated using **AliSim** (part of IQ-TREE 2) with parameters from Molloy and Warnow (2020).

Gene trees from SimPhy (Phase 2) are used directly - no branch length conversion is needed since SimPhy outputs branch lengths in **substitutions per site**.

### DNA Sequence Parameters

**Model**: GTR+Γ (General Time Reversible with Gamma rate heterogeneity)

**Substitution Rates** (AC, AG, AT, CG, CT, GT):
- Sampled from Dirichlet(12.776722, 20.869581, 5.647810, 9.863668, 30.679899, 3.199725)
- Independent sample for each gene tree

**Nucleotide Frequencies** (T, C, A, G):
- Sampled from Dirichlet(113.48869, 69.02545, 78.66144, 99.83793)
- Normalized to sum to 1.0
- Independent sample for each gene tree

**Gamma Shape Parameter** (α):
- Sampled from Lognormal(-0.470703916, 0.348667224)
- Controls rate heterogeneity across sites

**Alignment Length**: 1000 bp

### AliSim Command Example

```bash
# Sample parameters
gtr_rates="2.345678/3.456789/0.876543/1.234567/4.567890"  # Relative to GT=1.0
base_freqs="0.301822/0.211993/0.225259/0.260926"         # Sum to 1.0
alpha="0.879134"

# Build model and run
model="GTR{${gtr_rates}}+F{${base_freqs}}+G4{${alpha}}"

iqtree2 --alisim alignment \
  -m "${model}" \
  -t gene_tree.nwk \
  --length 1000 \
  -seed 22 \
  --out-format phylip
```

### Maximum Likelihood Gene Tree Estimation

Gene trees are estimated from true alignments using **PhyML** v3.1+ with:

**Model**: GTR+Γ
```bash
phyml -i alignment_TRUE.phy \
  -m GTR \
  -c 4 \
  -a e \
  -b 0 \
  -o tlr \
  --quiet
```

**Parameters:**
- `-m GTR`: General Time Reversible model
- `-c 4`: 4 discrete gamma rate categories
- `-a e`: Estimate gamma shape parameter
- `-b 0`: No bootstrap
- `-o tlr`: Optimize topology, branch lengths, and substitution rates

**Output**: Unrooted ML tree in `alignment_TRUE.phy_phyml_tree.txt`

**Note**: ML tree estimation can be skipped using the `-m` flag (alignments only), or run separately using the `-i` flag (infer-only mode). See "Two-Stage Workflow" section below for details.

### Expected ML Tree Quality

For 1000 bp alignments with appropriate parameters:
- **Tree size**: 5-15 substitutions per site
- **Parsimony score**: 1000-2500 substitutions
- **Log-likelihood**: Typically -5000 to -7000
- **Branch lengths**: Reasonable range (0.001-5.0)

## Phase 4: Simulating Protein Sequences and Estimating Gene Trees

Protein sequences are simulated similarly to DNA, but using amino acid substitution models.

### Protein Sequence Parameters

**Model**: WAG+Γ (Whelan and Goldman model with Gamma rate heterogeneity)

**Gamma Shape Parameter** (α):
- Sampled from Lognormal(-0.470703916, 0.348667224)
- Same distribution as DNA, independent samples

**Alignment Length**: 333 amino acids (≈ 1000 bp / 3)

### AliSim Command Example

```bash
# Sample gamma shape parameter
alpha="0.625"

# Build model and run
model="WAG+G4{${alpha}}"

iqtree2 --alisim alignment \
  -m "${model}" \
  -t gene_tree.nwk \
  --length 333 \
  -seed 22 \
  --out-format phylip
```

### Maximum Likelihood Estimation

Protein trees estimated using **PhyML** with:

```bash
phyml -i alignment_TRUE.phy \
  -d aa \
  -m WAG \
  -c 4 \
  -a e \
  -b 0 \
  -o tlr \
  --quiet
```

**Parameters:**
- `-d aa`: Amino acid sequences
- `-m WAG`: WAG substitution matrix
- Other parameters same as DNA

**Note**: ML tree estimation can be skipped using the `-m` flag (alignments only), or run separately using the `-i` flag (infer-only mode). See "Two-Stage Workflow" section below for details.

## Output Structure

For each replicate and configuration:

```
simulation_output/
├── species_trees/
│   └── lf12_replicate1.nex               # Species tree for Phase 2 (Nexus format)
├── gene_trees/
│   └── leaves12_dl1e-10_ps1e7/
│       └── replicate_1/
│           └── 1/
│               ├── s_tree.trees          # Species tree (from SimPhy)
│               ├── l_trees.trees         # Locus tree (after DL)
│               └── g_trees1.trees        # Gene tree (after ILS, branches in subs/site)
├── dna/
│   └── leaves12_dl1e-10_ps1e7/
│       └── replicate_1/
│           ├── gene_tree.nwk             # Gene tree from Phase 2 (branches in subs/site)
│           ├── alignment_TRUE.phy        # True alignment (from AliSim)
│           ├── alignment_TRUE.phy_phyml_tree.txt  # ML tree
│           ├── alignment_TRUE.phy_phyml_stats.txt # PhyML statistics
│           └── ml_gene_tree.nwk          # Copied ML tree
└── protein/
    └── leaves12_dl1e-10_ps1e7/
        └── replicate_1/
            ├── gene_tree.nwk             # Gene tree from Phase 2 (branches in subs/site)
            ├── alignment_TRUE.phy        # True alignment (from AliSim)
            ├── alignment_TRUE.phy_phyml_tree.txt  # ML tree
            ├── alignment_TRUE.phy_phyml_stats.txt # PhyML statistics
            └── ml_gene_tree.nwk          # Copied ML tree
```

## Simulation Scale

With default parameters:
- **Leaf-set sizes**: 2 (12 and 20 leaves)
- **DL rates**: 3 levels
- **Population sizes**: 2 levels
- **Replicates per configuration**: 100
- **Total configurations**: 2 × 3 × 2 = 12
- **Total datasets per sequence type**: 1,200
- **Total datasets (DNA + Protein)**: 2,400

## Pipeline Improvements (2025-12-23)

The pipeline was significantly improved to use SimPhy end-to-end:

### Previous Approach (TreeSim + SimPhy)
- Phase 1: TreeSim generated species trees
- Required scaling branch lengths (years → millions of years)
- Required scaling all rates by 10⁶
- Required regex manipulation to remove root branches
- Prone to scaling errors and tree format issues

### Current Approach (SimPhy Only)
- Phase 1: **SimPhy generates species trees directly**
- Phase 2: SimPhy generates gene trees (same tool, consistent behavior)
- No scaling needed - all parameters in years
- SimPhy outputs gene tree branch lengths directly in **substitutions per site**
- No tree manipulation required
- Simpler, more reliable, fewer sources of error

### Bug Fixes Applied
1. **Path mismatch**: Phase 2 and Phase 3 used inconsistent replicate directory names
2. **Incorrect scaling removal**: Removed obsolete branch length scaling (×0.0004) from Phase 3/4
3. **Octal number errors**: Fixed arithmetic with zero-padded numbers (01-09)

**⚠️ WARNING**: Any simulations generated with the old TreeSim-based pipeline should be regenerated.

See `BUGFIX_DL_AND_SEED.md` for documentation of earlier fixes.

## References

- Molloy, E.K. and Warnow, T. (2020). TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees.
- Rasmussen, M.D. and Kellis, M. (2011). A Bayesian approach for fast and accurate gene tree reconstruction.
- SimPhy: Mallo, D., De Oliveira Martins, L., and Posada, D. (2016). SimPhy: Phylogenomic simulation of gene, locus, and species trees. Systematic Biology, 65(2):334-344.
- **AliSim**: Ly-Trong, N. et al. (2022). AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the Genomic Era. Molecular Biology and Evolution, 39(5):msac092.
- **IQ-TREE 2**: Minh, B.Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5):1530-1534.
- PhyML: Guindon, S. et al. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology, 59(3):307-321.

## Command-Line Usage

The pipeline can be run with customizable parameters:

```bash
# Default: 100 replicates, 12 and 20 leaves, both DNA and protein
./simulate_pipeline.sh

# Custom configuration
./simulate_pipeline.sh -n 50 -l 12,20,50 -t both -o my_output -s 42

# Two-stage workflow: alignments first, ML trees later
# Stage 1: Generate alignments only (skip PhyML, ~3x faster)
./simulate_pipeline.sh -n 100 -l 12 -o my_output -m

# Stage 2: Infer ML trees from existing alignments
./simulate_pipeline.sh -o my_output -n 100 -l 12 -i

# Options:
#   -n NUM    Number of replicates per configuration (default: 100)
#   -l LEAVES Comma-separated leaf counts (default: 12,20)
#   -t TYPE   Sequence type: dna, protein, or both (default: both)
#   -o DIR    Output directory (default: simulation_output)
#   -s SEED   Random seed (default: 22)
#   -f NUM    Filter gene trees: require exactly NUM leaves (default: 0 = no restriction)
#   -m        Skip ML tree estimation (generate alignments only)
#   -i        Infer-only mode: run only PhyML on existing alignments
```

### Two-Stage Workflow

The pipeline supports a two-stage workflow for flexibility and efficiency:

**Stage 1: Simulation and Alignment Generation** (`-m` flag)
- Runs Phases 1-4 but **skips PhyML** tree inference
- Generates species trees, gene trees, and sequence alignments
- Approximately **3× faster** than full pipeline
- Useful for quickly generating alignment datasets

**Stage 2: ML Tree Inference** (`-i` flag)
- **Skips** Phases 1-2 (species and gene tree simulation)
- **Skips** sequence generation in Phases 3-4
- **Runs only PhyML** on existing `alignment_TRUE.phy` files
- Generates ML trees and PhyML statistics
- Must use **same parameters** (`-n`, `-l`) as Stage 1

**Use Cases:**
1. **Fast testing**: Generate alignments quickly, infer trees only if needed
2. **Batch processing**: Generate alignments once, distribute ML inference across compute nodes
3. **Parameter exploration**: Generate alignments once, experiment with different ML inference settings

### Gene Tree Filtering by Leaf Count

Due to duplication and loss events in Phase 2, gene trees can have varying numbers of leaves. The `-f NUM` option filters gene trees to ensure all replicates have exactly NUM leaves.

**Example:**
```bash
# Require all gene trees to have exactly 10 leaves
./simulate_pipeline.sh -n 100 -l 12 -f 10
```

**Filtering Process:**
1. SimPhy generates a gene tree with duplication/loss and ILS
2. The pipeline counts the number of leaves in the gene tree
3. If the tree has exactly NUM leaves, it's accepted for Phases 3-4
4. If not, SimPhy regenerates with a different seed (up to 1000 attempts)
5. Progress messages indicate how many retries were needed

**Technical Details:**
- Leaf counting uses a simple algorithm: count commas in Newick format + 1
- Each retry uses a different random seed: `seed + retry_count × 1,000,000`
- Maximum 1000 retries per replicate (configurable via `MAX_GENE_TREE_RETRIES`)
- If maximum retries reached, the last tree is accepted with a warning

**Example output:**
```
Configuration: leaves12_dl1e-10_ps1e7
  Replicate 1: accepted after 15 retries (10 leaves)
  Replicate 2: accepted after 3 retries (10 leaves)
  ...
```

**Use Cases:**
- **Consistent tree sizes**: Ensure all datasets have the same number of taxa
- **Fair comparison**: Remove variability due to different tree sizes
- **Controlled experiments**: Study effects of ILS and DL while maintaining constant tree size
- **Filtering artifacts**: Exclude gene trees with extreme duplication/loss outcomes

See `README.md` for complete usage documentation.
