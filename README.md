# HPVxHunter 

## Overview
HPVxHunter is an advanced bioinformatics pipeline engineered for Human Papillomavirus (HPV) sequence alignment, taxonomic classification, and novel lineage discovery. It implements an adaptive routing system to process curated reference databases and dynamic NCBI repositories, ensuring strict ICTV taxonomic compliance while capturing divergent sequences for downstream analysis.

## Core Capabilities

- **Adaptive Database Routing:**
  * **PaVE:** Executes parallel VSEARCH alignments against official ICTV reference types (`pave_classic`) and unofficial (`pave_non_official`). This segregation prevents misclassification, while keeping unofficial sequences available for valuable additional insights. Supports explicit database versioning (e.g., `--db pave_2026`).
  * **GenBank:** Standard alignment against broader NCBI archives. The `--update_genbank` flag forces an automatic API connection to download and replace the local archive with the latest sequences.
  * **Custom DB:** Accepts any user-provided FASTA file for alignment of unpublished sequences (allowing research teams to build local reference libraries while awaiting formal PaVE integration).
  
- **Precision Taxonomic Classification:** 
  * Direct lineage, sublineage, and high-risk flagging when using PaVE.
  * Two-step back-mapping for GenBank/Custom databases (identifies the best hit, then maps it against PaVE for robust annotation).
  * Built-in sequence length-based false-positive filtering.
  
- **Nomenclature Standardization:** HPVxHunter preserves official PaVE taxonomic designations. When a genotype has established lineages but a non-subdivided root clade in the reference header (e.g., `HPV16-A`), it is standardized to `A1` (`A-1`) to ensure consistent downstream plotting. However, for genotypes represented by a single reference genome without established lineages or sub-lineages (e.g., `HPV25REF`), the pipeline retains the `REF` designation (`Lineage: REF`, `Sub-lineage: REF`, or `Hit PAVE match: HPV25REF` in GenBank/Custom mode) rather than forcing an `A1` label, avoiding the misleading implication that other lineages or sub-lineages have been identified for that type.

- **Decoupled Novelty Detection (Deep Analysis):**
  * **Primary Screening & Thresholds:** Sequences exhibiting < 90.0% global whole-genome identity are assigned `ND` (Not Determined) in primary reports to prevent erroneous formal assignments, while still reporting the closest available database hit. Sequences with < 95.0% global identity are proactively flagged as divergent candidates.
  * **L1 ORF Verification:** When `--deep_analysis 1` is enabled, flagged sequences are routed to an integrated PuMA module that extracts the L1 open reading frame (ORF) and compares it against official PaVE L1 references. In strict compliance with ICTV criteria, candidate sequences displaying < 90.0% L1 nucleotide identity (> 10.0% divergence) are confirmed as novel HPV types in a dedicated secondary report without altering primary analysis reports.

---

## Output Interpretation & Key Columns

When running `HPVxHunter`, the output Excel report includes two critical columns designed to ensure robust quality control and taxonomic traceability across all database modes:

- **`Length Difference`** *(All modes)*: Measures the absolute base-pair length difference (in bp) between the query sequence and its best database reference hit. This serves as an essential Quality Control (QC) filter: large length deviations (e.g., >20-30 bp) indicate potential sequencing indels, homopolymeric errors, or incomplete genome assemblies that can artificially depress global similarity scores, preventing false-positive novel lineage or variant calls.
- **`Hit PAVE match`** *(GenBank and Custom Database modes only)*: When aligning against GenBank or custom FASTA files, the raw headers typically do not contain formal type, lineage, and sub-lineage annotations. To resolve this, `HPVxHunter` performs a two-step back-mapping: it identifies the best database hit, and then projects that reference sequence against the official PaVE database. The `Hit PAVE match` column displays this resulting official PaVE classification (Type-Lineage-Sublineage).
  * *How to read it:* While the `Similarity` and `Results` columns evaluate your query against its direct GenBank/custom hit, the `Hit PAVE match` column reveals what that hit officially represents in the PaVE nomenclature framework (e.g., if your query matches an unannotated GenBank entry that `Hit PAVE match` resolves as `HPV72-A-ND`, it confirms a real, circulating sub-lineage not yet formalized in PaVE).

### Example Output & Case Scenarios (GenBank / Custom DB Mode)

| Sequence | Match | Type | Lineage | Sub-lineage | Similarity | Length Difference | Results | Hit PAVE match |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `Sample_345` | `MF288727.1` | Same type | Same lineage | Same sublineage | 99.9% | 0 | Perfect or near-perfect match | `HPV18-A-A5` |
| `Sample_512` | `OP712107.1` | Same type | Same lineage | Same sublineage | 99.8% | 2 | Perfect or near-perfect match | `HPV72-A-ND` |
| `Sample_809` | `KC470285.1` | Same type | Same lineage | Different sublineage | 99.2% | 48 | Probably new sublineage | `HPV35-A-A1` |

- **Case 1 (`Sample_345` - Known sub-lineage):** The query matches GenBank entry `MF288727.1` with 99.9% identity and no length discrepancy (`0 bp`). Because `Hit PAVE match` confirms that `MF288727.1` corresponds to `HPV18-A-A5` in PaVE, `Sample_345` also belongs to the `HPV18-A-A5` sub-lineage.
- **Case 2 (`Sample_512` - Unassigned circulating sub-lineage):** The query is nearly identical (`99.8%`, `2 bp` difference) to GenBank isolate `OP712107.1`, meaning both sequences belong to the exact same sub-lineage. However, `Hit PAVE match` returns `HPV72-A-ND`, indicating that `OP712107.1` is not yet assigned to an official sub-lineage in PaVE. This tells the user that their sample belongs to a real, circulating HPV72 Lineage A variant already present in NCBI, even though it awaits formal PaVE sub-lineage assignment.
- **Case 3 (`Sample_809` - Sequencing / assembly artifact flagged by QC):** Based strictly on its `99.2%` similarity to `KC470285.1` (`Hit PAVE match`: `HPV35-A-A1`), the raw threshold flags `Sample_809` as a `"Probably new sublineage"`. However, the `Length Difference` column reveals a `48 bp` deviation, alerting the user that the depressed similarity score is likely driven by a technical deletion or incomplete genome assembly rather than true evolutionary divergence.
## Project Structure
```text
HPVxHunter/
├── core/
│   ├── analysis.py
│   ├── visualize.py
│   └── ...
├── data/                    
│   ├── pave_classic_2026.fasta
│   ├── pave_non_official_2026.fasta
│   ├── DB_full_L1.fasta
│   └── README.md
├── scripts/
│   └── puma/                # PuMA dependency for Deep Analysis (L1 extraction)
├── test_data/               # Provided sample data for testing
│   └── test.fasta
├── results/
├── main.py                # Main pipeline entry point
├── environment.yml          # Conda environment configuration
└── README.md
```

## Installation & Deployment

### Option A: Docker Container (Recommended)
Using the pre-compiled Docker image eliminates local dependency management (Conda, VSEARCH, Python). Docker Engine must be installed and running on the host system.

```bash
docker pull fmoa/hpvxhunter:latest
```

### Option B: Local Conda Environment
For active development or local execution without containerization.

```bash
git clone https://github.com/moa-f/HPVxHunter.git
cd HPVxHunter
conda env create -f environment.yml
conda activate HPVxHunter
```

## Usage via Docker

Due to container isolation, host directories must be explicitly mounted via the `-v` flag to allow I/O operations for input sequences and output results.

### 1. Environment Validation (Test Run)
The image natively embeds the `test_data` directory. To validate the container runtime and pipeline integrity, run the following command. Only the `results` volume requires mounting for data persistence:

```bash
docker run --rm \
  -v $(pwd)/results:/app/results \
  fmoa/hpvxhunter:latest \
  --db pave \
  --sample test_data/test.fasta \
  --run_name "validation_run"
```


### 2. Production Execution (Real Data)
Docker operates in an isolated file system. To allow the software to read your sequences, you must explicitly link a local directory to the container.

**Step 1:** Create a new directory named `data` in your current working location and place your FASTA file inside it (e.g., `clinical_01.fasta`).

**Step 2:** Run the pipeline using the following command. Modify `clinical_01.fasta` to match your actual file name, and adjust `clinical_cohort_A` for your desired output folder name.


```bash
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  fmoa/hpvxhunter:latest \
  --db genbank \
  --sample data/clinical_01.fasta \
  --run_name "clinical_cohort_A" \
  --deep_analysis 1
```

## Usage via Local Conda

If running directly on the host machine via the Conda environment, volume mounting is not required.

**1. Default Execution:**
```bash
python3 main.py --sample test_data/test.fasta --run_name "clinical_cohort_A"
```

**2. GenBank Execution with API Update:**
```bash
python3 main.py --db genbank --update_genbank --sample test_data/test.fasta
```

**3. Deep Analysis Execution:**
```bash
python3 main.py --db pave_2026 --sample test_data/test.fasta --deep_analysis 1
```
## Command-Line Arguments

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--db` | String | `pave` | Database selection (`pave`, `pave_XXXX`, `genbank`, or custom `.fasta` path). |
| `--sample` | String | `test.fasta` | Input sequences in FASTA format. |
| `--run_name` | String | Auto | Session identifier. Output files are saved in `analyze_results/<run_name>`. |
| `--update_genbank` | Flag | False | Forces NCBI API connection to download the latest GenBank sequences. |
| `--deep_analysis`| Integer | `0` | Set to `1` to enable L1 ORF extraction and novel type detection via PuMA. |
| `--plot` | String | `png` | Export format for final graphical distributions (`png` or `pdf`). |

## Acknowledgments & Citations
If you use HPVxHunter in your research, please ensure you also cite the foundational tools it relies upon:
* **VSEARCH:** Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584.
* **PaVE (Papillomavirus Episteme):** Van Doorslaer K, et al. (2017). The Papillomavirus Episteme: a major update to the papillomavirus sequence database. *Nucleic Acids Res*.
* **PuMA:** Pace J, Youens-Clark K, Freeman C, Van Doorslaer K. (2020). PuMA: A papillomavirus genome annotation tool. *Virus Evolution*, 6(2), veaa068.
