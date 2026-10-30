# MitoArrange

**Comparative analysis and visualization of mitochondrial gene organization from annotated genomes**

MitoArrange compares **gene content, circular gene order, and gene orientation** across mitochondrial genomes using GenBank annotations. It identifies conserved gene adjacencies, differences in gene arrangements, conserved blocks, and gene-order patterns compatible with a single inversion. It also groups genomes by identical arrangements or by an adjustable similarity threshold and generates publication-oriented SVG diagrams.

**MitoArrange is an annotation-based gene-arrangement comparison tool, not a nucleotide-level synteny aligner or a read-supported structural-variant caller.** It does not reconstruct physical mitochondrial DNA molecules.

## Features

- Reads `.gb`, `.gbk`, `.gbf`, and `.genbank` files from a directory (recursively), a ZIP archive, or an individual file.
- Standardizes the displayed circular order using a start/orientation gene (default: `cox1`).
- Compares **signed circular gene adjacencies** (gene order plus strand orientation), independent of circular rotation and whole-sequence reverse complementation.
- Reports conserved blocks, changed gene junctions, pairwise breakpoint distances, and candidate single-block inversions.
- Filters repeated genes, discontinuous/multipart genes, and genes with unknown strand; recognizes contiguous genes crossing a circular sequence origin.
- Allows a minimum gene-presence frequency (`--min_gene_presence`, default: `1.0`), while keeping the other exclusion rules strict.
- Groups samples by input folder, exact gene arrangement, or complete-linkage similarity clustering.
- Exports detailed tab-separated tables (`.tsv`) and standard/compact SVG figures with stable gene colors and textures.
- Requires **no sequence alignment, read mapping, third-party command-line executables, or SciPy**.

## Installation

1. Install **Python 3.10 or newer**.
2. Install Biopython:

   ```bash
   python -m pip install biopython
   ```

   Alternatively, if you have cloned the repository:

   ```bash
   python -m pip install -r requirements.txt
   ```

3. Download `mitoarrange.py` from this repository and run:

   ```bash
   python mitoarrange.py --help
   ```

> On some systems, use `python3` instead of `python`.

## Quick start

Compare all annotated genomes in a folder:

```bash
python mitoarrange.py -i genomes/ -o results/
```

Select a reference and an orientation gene:

```bash
python mitoarrange.py \
  -i genomes/ \
  -o results/ \
  --reference Bacuri_Acu \
  --start-gene cox1
```

Read a ZIP archive, group genomes by similarity, and show every member's name beside the representative:

```bash
python mitoarrange.py \
  -i genomes.zip \
  -o results/ \
  --group-by similarity \
  --max-breakpoints 3 \
  --collapsed_sub_all
```

Retain genes annotated in at least 80% of samples, using larger text and wider arrows:

```bash
python mitoarrange.py \
  -i genomes/ \
  -o results/ \
  --min_gene_presence 0.8 \
  --text_scale 1.3 \
  --arrow_width_scale 1.2 \
  --arrow_height_scale 1.2
```

If `-o` is omitted, the output folder is `mitoarrange_results/`.


## How the analysis works

### 1. Select comparable genes

For each gene name, MitoArrange evaluates its presence across the entire dataset and applies global exclusions:

| Condition | Behavior |
|---|---|
| Missing from some genomes | Retained only if its fraction of present samples meets `--min_gene_presence` (default `1.0`). |
| More than one `gene` feature with the same name in any sample | Globally excluded as `multicopy`. |
| Discontinuous/multipart `gene` feature in any sample | Globally excluded as `multipart`, except the specific circular-origin case below. |
| Unknown/undefined gene strand in any sample | Globally excluded as `unknown_strand`. |
| Two-part join spanning positions at the end and beginning of a circular record | Treated as one contiguous gene, **not** automatically as a trans-spliced gene. |

The circular-origin exception applies only when a gene has **exactly two parts**, one touching the sequence end and the other touching position zero. Genes with three or more parts are not granted that exception.

For example, `--min_gene_presence 0.8` retains a valid gene present in at least 8 of 10 samples; it **does not** override a multicopy, multipart, or unknown-strand exclusion. The `gene_presence.tsv` and `excluded_genes.tsv` files document these decisions.

### 2. Normalize the displayed gene order

The start gene (`cox1` by default) defines the displayed orientation and starting point. If it is available as one gene with known strand, it may be used for orientation even when excluded from the comparative marker set. Otherwise, MitoArrange falls back to the alphabetically first retained gene in that sample.

This normalization helps make diagrams comparable; the underlying signed-adjacency representation is independent of where a circular genome sequence starts or whether the whole sequence is reverse-complemented.

### 3. Compare signed circular adjacencies

Each retained gene has a name and strand (`+` or `-`). MitoArrange constructs the neighboring **head/tail extremity adjacencies** around the circular order and compares those sets between genomes.

- `changed_adjacencies` counts adjacencies present in the first genome but not the second (after any required shared-gene projection).
- `normalized_distance` divides this count by the number of genes shared by that genome pair.
- Conserved blocks are uninterrupted runs of shared adjacencies relative to the reference order.
- A candidate **single exact signed inversion** is reported only when reversing one contiguous signed-gene interval reproduces the entire compared order and there are exactly two changed adjacencies.

When the minimum presence threshold is below `1.0`, each pairwise comparison is projected onto genes found in **both** samples. This avoids treating a missing annotation as a breakpoint by itself. However, reduced gene overlap also reduces structural resolution.

**Important:** A changed adjacency is not necessarily one DNA breakpoint or one evolutionary rearrangement event. The script cannot determine exact nucleotide breakpoints from gene annotations.

### 4. Organize samples for visualization

Use `--group-by` to choose how visualization groups are defined:

| Mode | Meaning |
|---|---|
| `auto` (default) | Use source folders when at least two distinct source-folder groups occur; otherwise use exact gene-arrangement groups. |
| `folder` | Group by the input file's immediate parent folder grouping. Use this if you decide to manually sort the groups, create a folder for each group and put your annotations files in those folders. |
| `structure` | Group samples with identical complete signed circular adjacency signatures (`S01`, `S02`, ...). |
| `similarity` | Group by pairwise breakpoint distance using deterministic **complete-linkage** agglomerative clustering (`C01`, `C02`, ...). |

Similarity grouping requires explicitly choosing `--group-by similarity`. By default its maximum allowed distance is zero; provide **one** of the following cutoffs:

```bash
--group-by similarity --max-breakpoints 3
```

or

```bash
--group-by similarity --max-breakpoint-fraction 0.15
```

With complete linkage, all pairwise distances within a completed group must be no greater than the specified cutoff. The fraction uses each pair's number of shared retained genes as its denominator.

**Exact structure groups** require identical full gene/adjoining-extremity signatures. If a gene is missing from one genome, two samples can be in different strict `S` groups even when their shared-gene-projected pairwise distance is zero.

### 5. Choose representatives and row order

`--representative` chooses the actual genome displayed for each collapsed group:

- `auto` (default): **medoid** for similarity groups; alphabetical first sample for other modes.
- `medoid`: sample with the smallest summed within-group breakpoint distance (alphabetical tie-breaker).
- `alphabetical`: alphabetically first sample.


`--display_order similarity` (default) arranges groups and sample rows heuristically to put similar gene orders near each other. `--display_order alphabetical` switches to alphabetic ordering. **Row ordering is not phylogenetic inference** and does not change group membership.

## Visualization

MitoArrange generates two SVG styles:

- **Full gene-arrow view:** an arrow per gene, labeled with its name; arrow direction represents relative strand orientation.
- **Compact view:** consecutive same-strand genes are combined into directional blocks; the internal colored/hatched bars distinguish individual genes. The compact view is intended to be included in the article main text.

Both views use a shared input-wide gene identity palette, beginning with 12 contrasting colors and adding textures for further genes. SVGs are editable in vector graphics software. Arrow size represents the schematic gene arrangement, **not nucleotide length**.

You can adjust figure dimensions:

| Option | Default | Effect |
|---|---:|---|
| `--text_scale` | `1.0` | Scales figure typography and associated label/legend spacing. |
| `--arrow_width_scale` | `1.0` | Scales gene-arrow width and compact tracks. |
| `--arrow_height_scale` | `1.0` | Scales arrow height and the compact view's internal color-bar height. |
| `--collapsed_sub_all` | off | Shows all members' names alongside the representative in collapsed diagrams. |

For large groups, `--collapsed_sub_all` can create long labels; you may need to edit their placement in your final figure.

## Output files

The output directory contains:

```text
results/
├── README_results.md
├── overview.svg
├── overview_compact_view.svg
├── overview_collapsed.svg
├── overview_collapsed_compact_view.svg
├── individual_image_groups/
│   ├── <group>.svg
│   └── <group>_compact_view.svg
├── samples.tsv
├── retained_genes.tsv
├── gene_order.tsv
├── gene_presence.tsv
├── excluded_genes.tsv
├── reference_comparisons.tsv
├── pairwise_comparisons.tsv
├── pairwise_blocks.tsv
├── changed_junctions.tsv
├── exact_inversions.tsv
├── all_pairs.tsv
├── breakpoint_distance_matrix.tsv
├── breakpoint_distance_pairs.tsv
├── visualization_groups.tsv
├── cluster_membership.tsv
├── cluster_summary.tsv
└── display_order.tsv
```

SVG files inside `individual_image_groups/` are generated for each visualization group, with corresponding compact versions. The actual group filenames depend on the chosen grouping and sample metadata.

| File | Contents |
|---|---|
| `README_results.md` | Run-specific methods summary, group memberships, reference comparisons, interpretation notes, and file inventory. |
| `samples.tsv` | Input genome metadata, annotation counts, retained-gene counts, anchor information, and strict structure IDs. |
| `retained_genes.tsv` | Genes selected for analysis across the dataset. |
| `gene_order.tsv` | Normalized retained gene order, strands, and source coordinates for each sample. |
| `gene_presence.tsv` | Presence frequency, retention decision, missing samples, and exclusion reasons for each gene. |
| `excluded_genes.tsv` | Genes rejected from comparisons and the reasons for exclusion. |
| `reference_comparisons.tsv` | Gene-adjacency comparison against the selected reference. |
| `pairwise_comparisons.tsv` | One comparison per unordered genome pair, including conserved blocks and classifications. |
| `pairwise_blocks.tsv` | Conserved-gene blocks, positions, and relative orientations. |
| `changed_junctions.tsv` | Head/tail gene junctions unique to either member of each compared pair. |
| `exact_inversions.tsv` | Intervals compatible with one exact signed-gene inversion. |
| `all_pairs.tsv` | All ordered genome pairs including self comparisons, strict structure agreement, and shared-gene counts. |
| `breakpoint_distance_matrix.tsv` | Square matrix of pairwise changed-adjacency counts. |
| `breakpoint_distance_pairs.tsv` | Long-format pairwise distances and normalized distances. |
| `visualization_groups.tsv` / `cluster_membership.tsv` | Visualization membership, representatives, group statistics, and associated SVG paths. |
| `cluster_summary.tsv` | Per-group size, representative selection, and maximum within-group distance. |
| `display_order.tsv` | Sample and group row order used in exported SVGs. |

## Command-line reference

| Argument | Description |
|---|---|
| `-i, --input` | **Required.** GenBank file, recursively scanned folder, or ZIP of GenBank files. |
| `-o, --output` | Output directory (default: `mitoarrange_results`). |
| `--reference` | Reference sample filename stem or full filename (default: first alphabetic sample). |
| `--start-gene` | Gene for orientation/start (default: `cox1`). |
| `--min_gene_presence`, `--min-gene-presence` | Valid fraction greater than zero and at most one (default: `1.0`). |
| `--group-by` | `auto`, `folder`, `structure`, or `similarity` (default: `auto`). |
| `--max-breakpoints` | Nonnegative integer maximum changed adjacencies for similarity grouping (default: `0` if applicable). |
| `--max-breakpoint-fraction` | Fraction cutoff from `0` to `1` for similarity grouping; mutually exclusive with `--max-breakpoints`. |
| `--representative` | `auto`, `alphabetical`, or `medoid` (default: `auto`). |
| `--collapsed_sub_all` | Label collapsed representations with all group member names. |
| `--display_order`, `--display-order` | `similarity` (default) or `alphabetical`. |
| `--text_scale` | Positive text scaling factor (default: `1.0`). |
| `--arrow_width_scale` | Positive arrow-width scaling factor (default: `1.0`). |
| `--arrow_height_scale` | Positive arrow-height scaling factor (default: `1.0`). |

Run `python mitoarrange.py --help` for the built-in command help.

## Interpretation and limitations

MitoArrange compares **annotated arrangements of homologously named genes**, not aligned genomic sequences. Its results are best treated as hypotheses for further investigation, especially in plant mitochondria where repeats, recombination, alternative genomic configurations, and assembly choices may affect apparent arrangements.

- The algorithm does not validate rearrangements with long reads, read mapping, genome alignments, or independent assemblies.
- It cannot locate nucleotide-level rearrangement breakpoints. A changed junction may span a large intergenic region.
- Rearrangements that do not affect the analyzed genes or their adjacencies may be invisible.
- Missing genes and conservative exclusions can reduce resolution; a zero distance on a shared subset does not prove globally identical mitochondrial organization.
- A conserved block with reversed orientation is **not automatically** evidence for a single historical inversion.
- Comparisons depend on annotation correctness and consistent gene naming; MitoArrange is not an orthology inference tool.
- The diagrams show linearized representations of circular gene orders, **not physically validated mitochondrial conformations**.

## Citation

The paper that first applies this tool in under development.
