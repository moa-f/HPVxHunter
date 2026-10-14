# HPVxHunter 

## Overview
HPVxHunter is an advanced bioinformatics pipeline engineered for Human Papillomavirus (HPV) sequence alignment, taxonomic classification, and novel lineage discovery. It implements an adaptive routing system to process curated reference databases and dynamic NCBI repositories, ensuring strict ICTV taxonomic compliance while capturing divergent sequences for downstream phylogenetic mapping.

## Core Capabilities

* **Adaptive Database Routing:**
  * **PaVE:** Executes parallel VSEARCH alignments against official ICTV reference types (`pave_classic`) and unofficial  (`pave_non_official`). This segregation prevents misclassification. Supports explicit database versioning (e.g., `--db pave_2026`).
  * **GenBank:** Standard alignment against broader NCBI archives. The `--update_genbank` flag forces an automatic API connection to download and replace the local archive with the latest sequences.
  * **Custom DB:** Accepts any user-provided FASTA file for alignment of unpublished sequences.
  
* **Precision Taxonomic Classification:** 
  * Direct lineage, sublineage, and high-risk flagging when using PaVE.
  * Two-step back-mapping for GenBank/Custom databases (identifies the best hit, then maps it against PaVE for robust annotation).
  * Built-in sequence length-based false-positive filtering.
  
* **Nomenclature Standardization:** 
  * Isolated root lineages (e.g., Lineage "A") are automatically standardized to "A1" internally. This ensures absolute consistency across Excel reports and prevents fragmented sub-lineage distribution plots.
  
* **Decoupled Novelty Detection (Deep Analysis):**
  * **Strict ICTV Reporting:** Sequences with < 90.0% global identity are assigned `ND` (Not Determined) in primary reports, strictly adhering to the official definition of novel HPV types.
  * **Divergent Lineage Capture:** Sequences with < 95.0% identity silently trigger the `--deep_analysis` module. This extracts the L1 ORF using PuMA and compares it against PaVE's L1 references to detect emerging or moderately divergent lineages without corrupting the main Excel report.

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
├── main_2.py                # Main pipeline entry point
├── environment.yml          # Conda environment configuration
└── README.md
```

## Installation

We recommend using Conda to manage the required dependencies (Python, Pandas, Biopython, VSEARCH, etc.).

```bash
# Clone the repository
git clone https://github.com/moa-f/HPVxHunter.git
cd HPVxHunter

# Create and activate the conda environment
conda env create -f environment.yml
conda activate hpvxhunter
```
## Usage Examples

**1. Default Execution (PaVE 2026 Basic-Search):**
```bash
python3 main_2.py --sample test_data/test.fasta --run_name "clinical_cohort_A"
```

**2. GenBank Execution with API Update:**
```bash
python3 main_2.py --db genbank --update_genbank --sample test_data/test.fasta
```
**2. Deep Analysis Execution :**
```bash
python3 main_2.py --db pave_2026 --sample test_data/test.fasta --deep_analysis 1
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
