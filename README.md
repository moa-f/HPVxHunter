# HPVxHunter 

## Overview
HPVxHunter is an advanced bioinformatics pipeline engineered for Human Papillomavirus (HPV) sequence alignment, taxonomic classification, and novel lineage discovery. It implements an adaptive routing system to process curated reference databases and dynamic NCBI repositories, ensuring strict ICTV taxonomic compliance while capturing divergent sequences for downstream analysis.

## Core Capabilities

* **Adaptive Database Routing:**
  * **PaVE:** Executes parallel VSEARCH alignments against official ICTV reference types (`pave_classic`) and unofficial  (`pave_non_official`). This segregation prevents misclassification, while keeping unofficial sequences available for valuable additional insights. Supports explicit database versioning (e.g., `--db pave_2026`).
  * **GenBank:** Standard alignment against broader NCBI archives. The `--update_genbank` flag forces an automatic API connection to download and replace the local archive with the latest sequences.
  * **Custom DB:** Accepts any user-provided FASTA file for alignment of unpublished sequences.
  
* **Precision Taxonomic Classification:** 
  * Direct lineage, sublineage, and high-risk flagging when using PaVE.
  * Two-step back-mapping for GenBank/Custom databases (identifies the best hit, then maps it against PaVE for robust annotation).
  * Built-in sequence length-based false-positive filtering.
  
* **Nomenclature Standardization:** Isolated root lineages (e.g., Lineage "A") are automatically standardized to "A1". This ensures absolute consistency across Excel reports and prevents fragmented sub-lineage distribution plots.
* **Decoupled Novelty Detection (Deep Analysis):**
  * **Strict ICTV Reporting:** Sequences with < 90.0% global identity are assigned `ND` (Not Determined) in primary reports, strictly adhering to the official definition of novel HPV types.
  * **Divergent Lineage Capture:** Sequences with < 95.0% identity (an arbitrary threshold designed to proactively flag potential divergence within the PaVE nomenclature framework)  are marked for further investigation. When the `--deep_analysis 1` is enabled, the pipeline extracts the L1 ORF using PuMA and compares it against PaVE's L1 references to detect emerging lineages without altering the main Excel report.    
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
  --deep_analysis
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
