# Command-Line Options Reference

Both `simulate_pipeline.sh` (Bash) and `simulate_pipeline.py` (Python) now support flexible command-line options for customizing your simulation runs.

## Quick Reference

| Parameter | Bash | Python | Description | Default |
|-----------|------|--------|-------------|---------|
| **Replicates** | `-n NUM` | `-n, --replicates NUM` | Number of replicates per configuration | 100 |
| **Leaves** | `-l LEAVES` | `-l, --leaves LEAVES` | Leaf counts (comma-sep for bash, space-sep for python) | 12,20 or 12 20 |
| **Type** | `-t TYPE` | `-t, --type TYPE` | Sequence type: `dna`, `protein`, or `both` | `both` |
| **Output** | `-o DIR` | `-o, --output DIR` | Output directory path | `simulation_output` |
| **Seed** | `-s SEED` | `-s, --seed SEED` | Random seed for reproducibility | 22 |
| **ML Mode** | `-a MODE` | N/A | ML tree inference mode: `true`, `inferred`, `true,inferred`, `none` | `inferred` |
| **DL Model** | `--dl-model` | N/A | Duplication/loss rates: `low`, `medium`, `high`, or list | `low,medium,high` |
| **ILS Model** | `--ils-model` | N/A | ILS/population: `low`, `medium`, or list | `low,medium` |
| **Help** | `-h` | `-h, --help` | Display help message | - |

## Detailed Options

### -n, --replicates NUM

Number of replicate datasets to generate per configuration.

**Default:** 100

**Examples:**
```bash
# Quick test with 10 replicates
./simulate_pipeline.sh -n 10
./simulate_pipeline.py -n 10

# Large study with 200 replicates
./simulate_pipeline.sh -n 200
./simulate_pipeline.py --replicates 200
```

**Impact on output:**
- Each configuration generates NUM datasets
- Total datasets = NUM × (leaf_sizes × DL_rates × pop_sizes)
- Default: 100 × (2 × 3 × 2) = 1,200 datasets per sequence type

### -l, --leaves LEAVES

Species tree leaf counts to simulate.

**Default:** 12,20 (bash) or 12 20 (python)

**Bash format:** Comma-separated list
```bash
./simulate_pipeline.sh -l 12
./simulate_pipeline.sh -l 12,20
./simulate_pipeline.sh -l 12,20,50
```

**Python format:** Space-separated list
```bash
./simulate_pipeline.py -l 12
./simulate_pipeline.py -l 12 20
./simulate_pipeline.py --leaves 12 20 50
```

**Common values:**
- **12** - Small phylogeny (fast)
- **20** - Medium phylogeny (default)
- **50** - Large phylogeny (slower)
- **100** - Very large phylogeny (computationally intensive)

**Impact on output:**
- More leaves = larger trees
- Computational time scales approximately O(n²) with tree size
- More leaves = more realistic but slower simulation

### -t, --type TYPE

Which sequence types to simulate.

**Default:** both

**Options:**
- `dna` - DNA sequences only (Phase 3 only)
- `protein` - Protein sequences only (Phase 4 only)
- `both` - Both DNA and protein sequences (Phases 3 & 4)

**Examples:**
```bash
# DNA only (faster, ~50% time savings)
./simulate_pipeline.sh -t dna
./simulate_pipeline.py --type dna

# Protein only
./simulate_pipeline.sh -t protein
./simulate_pipeline.py -t protein

# Both (default)
./simulate_pipeline.sh -t both
./simulate_pipeline.py --type both
```

**Use cases:**
- Use `-t dna` when only testing DNA-based phylogenetic methods
- Use `-t protein` when only interested in protein evolution
- Use `-t both` for comprehensive comparison studies

**Runtime impact:**
- `dna` only: ~50% faster than `both`
- `protein` only: ~50% faster than `both`
- `both`: Full runtime but generates both data types

### -o, --output DIR

Output directory for all generated files.

**Default:** simulation_output

**Examples:**
```bash
# Custom output directory
./simulate_pipeline.sh -o my_simulation
./simulate_pipeline.py --output my_simulation

# Organized by date
./simulate_pipeline.sh -o results_2024_03_15
./simulate_pipeline.py -o results/experiment1
```

**Directory structure created:**
```
<output_dir>/
├── species_trees/
├── gene_trees/
├── dna/           # Only if -t dna or -t both
└── protein/       # Only if -t protein or -t both
```

### -s, --seed SEED

Random seed for reproducibility.

**Default:** 22

**Examples:**
```bash
# Use default seed (reproducible)
./simulate_pipeline.sh

# Custom seed
./simulate_pipeline.sh -s 42
./simulate_pipeline.py --seed 12345

# Different runs with different seeds
./simulate_pipeline.sh -s 1 -o run1
./simulate_pipeline.sh -s 2 -o run2
./simulate_pipeline.sh -s 3 -o run3
```

**Use cases:**
- Same seed = identical results (for reproducibility)
- Different seeds = independent replicate runs
- Useful for testing method variance across different random datasets

### -a MODE

ML tree inference mode for controlling which alignments to use for PhyML tree estimation.

**Default:** inferred

**Options:**
- `true` - Estimate ML trees only from true alignments (no MAFFT alignment inference)
- `inferred` - Estimate ML trees only from MAFFT-inferred alignments (default)
- `true,inferred` - Estimate ML trees from both true and inferred alignments (creates two sets of trees)
- `none` - Skip all ML tree estimation (same as `-m` flag)

**Examples:**
```bash
# Default: only inferred alignments
./simulate_pipeline.sh -n 10 -l 12 -a inferred

# Only true alignments, no MAFFT
./simulate_pipeline.sh -n 10 -l 12 -a true

# Both alignment types for comparison
./simulate_pipeline.sh -n 10 -l 12 -a true,inferred

# No ML trees at all
./simulate_pipeline.sh -n 10 -l 12 -a none
```

**File Output:**
- `-a true`: Creates `ml_gene_tree.nwk`
- `-a inferred`: Creates `ml_gene_tree.nwk` and `alignment_INFERRED.phy`
- `-a true,inferred`: Creates `ml_gene_tree_TRUE.nwk` and `ml_gene_tree_INFERRED.nwk`
- `-a none`: No ML tree files created

**Use Cases:**
- Use `inferred` to test phylogenetic methods on realistic (noisy) alignments
- Use `true` for perfect alignment scenarios or when MAFFT is not needed
- Use `true,inferred` to compare method accuracy between true and inferred alignments
- Use `none` when you only need simulated sequences without tree inference

### --dl-model MODEL

Duplication/loss rate preset for controlling gene family evolution.

**Default:** low,medium,high (all three rates)

**Options:**
- `low` - 1e-10 events/year (minimal gene family evolution)
- `medium` - 2e-10 events/year (moderate duplication/loss)
- `high` - 5e-10 events/year (high duplication/loss rate)
- Comma-separated list - e.g., `low,high` for only those two rates

**Examples:**
```bash
# Use only high duplication/loss rate
./simulate_pipeline.sh --dl-model high

# Use only low and medium rates
./simulate_pipeline.sh --dl-model low,medium

# Default: all three rates
./simulate_pipeline.sh --dl-model low,medium,high
./simulate_pipeline.sh  # Same as above
```

**Impact on configurations:**
- Each DL rate creates a separate configuration
- Total configs = #leaves × #DL_rates × #ILS_levels
- Example: `-l 12 --dl-model high --ils-model low` = 1 configuration
- Example: `-l 12,20 --dl-model low,high --ils-model low,medium` = 8 configurations

**Use Cases:**
- Use `--dl-model high` to focus on challenging gene family scenarios
- Use `--dl-model low` for minimal duplication/loss (closer to species tree)
- Use multiple rates to study method robustness across evolutionary scenarios

### --ils-model MODEL

ILS (Incomplete Lineage Sorting) model controlling coalescent variation via population size.

**Default:** low,medium (both levels)

**Options:**
- `low` - Population size 1e7 (high ILS, smaller population, more coalescent variation)
- `medium` - Population size 5e7 (low ILS, larger population, less coalescent variation)
- Comma-separated list - e.g., `low,medium` for both

**Examples:**
```bash
# Use only high ILS (smaller population)
./simulate_pipeline.sh --ils-model low

# Use only low ILS (larger population)
./simulate_pipeline.sh --ils-model medium

# Default: both ILS levels
./simulate_pipeline.sh --ils-model low,medium
./simulate_pipeline.sh  # Same as above
```

**Biological Meaning:**
- **Low ILS (`medium` population)**: Larger effective population → more time for lineages to coalesce → gene trees more similar to species tree
- **High ILS (`low` population)**: Smaller effective population → less time for coalescent → gene trees more discordant from species tree

**Impact on configurations:**
- Each ILS level creates a separate configuration
- Total configs = #leaves × #DL_rates × #ILS_levels
- Example: `--dl-model low --ils-model medium` = 1 configuration per leaf count
- Example: `--dl-model low,high --ils-model low,medium` = 4 configurations per leaf count

**Use Cases:**
- Use `--ils-model low` to study high coalescent variation (challenging for species tree methods)
- Use `--ils-model medium` for more concordant gene-species tree relationships
- Use both to evaluate method performance across ILS levels

## Common Usage Patterns

### Quick Test Run
```bash
# Fast test: 10 replicates, 12 leaves, DNA only
./simulate_pipeline.sh -n 10 -l 12 -t dna
./simulate_pipeline.py -n 10 -l 12 -t dna
```
**Time:** ~5-10 minutes
**Output:** 60 DNA datasets (10 × 6 configurations)

### Small Study
```bash
# 50 replicates, default leaves, both types
./simulate_pipeline.sh -n 50
./simulate_pipeline.py -n 50
```
**Time:** ~2-4 hours
**Output:** 600 DNA + 600 protein datasets

### Default (Full Study)
```bash
# 100 replicates, 12 and 20 leaves, both types
./simulate_pipeline.sh
./simulate_pipeline.py
```
**Time:** ~5-9 hours
**Output:** 1,200 DNA + 1,200 protein datasets

### Large Phylogeny Study
```bash
# 100 replicates, 50 leaves, both types
./simulate_pipeline.sh -l 50
./simulate_pipeline.py -l 50
```
**Time:** ~12-24 hours
**Output:** 600 DNA + 600 protein datasets with 50-taxon trees

### Multi-size Comparison
```bash
# Compare small, medium, and large trees
# Bash
./simulate_pipeline.sh -l 12,20,50 -n 50

# Python
./simulate_pipeline.py -l 12 20 50 -n 50
```
**Time:** ~8-16 hours
**Output:** 900 DNA + 900 protein datasets across 3 tree sizes

### DNA vs Protein Performance Test
```bash
# Generate DNA data
./simulate_pipeline.sh -t dna -o dna_only

# Generate protein data (different run)
./simulate_pipeline.sh -t protein -o protein_only
```
**Time:** ~2.5-4.5 hours each
**Output:** Separate datasets for comparing DNA vs protein methods

### Reproducibility Test
```bash
# Run 1 with seed 22
./simulate_pipeline.sh -s 22 -o run1

# Run 2 with same seed (should be identical)
./simulate_pipeline.sh -s 22 -o run2

# Verify identical output
diff -r run1/ run2/  # Should show no differences
```

### Independent Replicates
```bash
# Generate 3 independent runs with different seeds
./simulate_pipeline.sh -s 1 -o replicate1 &
./simulate_pipeline.sh -s 2 -o replicate2 &
./simulate_pipeline.sh -s 3 -o replicate3 &
wait
```

## Configuration Matrix

The total number of configurations depends on:
- Number of leaf counts specified
- Fixed DL rates: 3 (10⁻¹⁰, 2×10⁻¹⁰, 5×10⁻¹⁰)
- Fixed population sizes: 2 (10⁷, 5×10⁷)

**Formula:** Total configs = `#leaves × 3 × 2`

**Examples:**

| Leaves | Configs | Total Datasets (per type) |
|--------|---------|---------------------------|
| -l 12 | 6 | 6 × NUM_REPLICATES |
| -l 12,20 | 12 | 12 × NUM_REPLICATES |
| -l 12,20,50 | 18 | 18 × NUM_REPLICATES |

Default (`-l 12,20 -n 100`): **1,200 datasets per sequence type**

## Tips and Best Practices

1. **Start small:** Always test with `-n 10 -l 12 -t dna` first
2. **Use specific types:** If you only need DNA, use `-t dna` to save time
3. **Organize output:** Use `-o` to create descriptive output directories
4. **Multiple runs:** Use different seeds (`-s`) for independent replicates
5. **Reproducibility:** Document the exact command and seed used
6. **Disk space:** Check available space before large runs (use `df -h`)
7. **Monitoring:** Use `watch -n 60 'ls -lh simulation_output/*/*/replicate*'` to monitor progress

## Troubleshooting

**"Invalid sequence type" error:**
```bash
# Wrong
./simulate_pipeline.sh -t DNA  # Case matters!

# Correct
./simulate_pipeline.sh -t dna
```

**Bash leaf separator:**
```bash
# Wrong (bash)
./simulate_pipeline.sh -l 12 20  # Space doesn't work in bash

# Correct (bash)
./simulate_pipeline.sh -l 12,20  # Use comma
```

**Python leaf separator:**
```bash
# Wrong (python)
./simulate_pipeline.py -l 12,20  # Comma doesn't work in python

# Correct (python)
./simulate_pipeline.py -l 12 20  # Use space
```

## See Also

- `README.md` - Full documentation
- `QUICK_START.md` - Quick start guide
- `config.example.sh` - Advanced configuration template
- `DNA_vs_PROTEIN.md` - Comparison of sequence types
