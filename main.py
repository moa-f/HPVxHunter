import argparse
import os
import re
import sys
import subprocess
import datetime
from pathlib import Path
import pandas as pd 
from Bio import Entrez
import time
import logging
from core.analysis import (process_vsearch_clusters, 
                           generate_final_custom_report, 
                           extract_fasta_sequences, 
                           prepare_deep_analysis_fasta,
                           extract_l1_sequences,
                           finalize_deep_analysis,
                           sanitize_sample_fasta,
                           DataValidationError)

from core.visualize import generate_plots

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("EvolutionaryPipeline")

EXCLUSION_PATTERN = re.compile(
    r'patent|-A/|vaccine|method|quantification|amplification|diagnostic|transcript|splice|primer|probe',
    re.IGNORECASE
)

SCRIPT_DIR = Path(__file__).parent.resolve()

def parse_arguments():
    # Détermination du chemin absolu du répertoire contenant le script 

    parser = argparse.ArgumentParser(
        description="Computational pipeline for discovering novel human papillomavirus (HPV) sub-lineages.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Construction des chemins par défaut absolus
    default_db_path = SCRIPT_DIR / "data" / "full_database_Pave.fasta"
    default_sample_path = SCRIPT_DIR / "test_data" / "test.fasta"

    # Paramètre : --db
    parser.add_argument("--db", type=str, default=default_db_path,
                        help="Reference database to use for the analysis.\n"
                             " - 'pave': loads the default version (2026).\n"
                             " - 'pave_XXXX': loads a specific version (e.g., pave_2026).\n"
                             " - 'genbank': uses the latest downloaded NCBI archive.\n"
                             " - Or specify an absolute path to a custom FASTA file.")
    
    
    # Paramètre : ----update_genbank
    parser.add_argument("--update_genbank", action="store_true",
                        help="Forces the pipeline to connect to NCBI API and download the latest sequences into a new timestamped archive.")

    # Paramètre : --sample
    parser.add_argument("--sample", type=str, default=default_sample_path,
                        help="Input sample sequences in FASTA format.\n"
                             "Default: test_data/test.fasta")

    # Paramètre : --run_name 
    parser.add_argument("--run_name", type=str, default=None,
                        help="Name of the run/session. Used for naming outputs.\n"
                             "If not specified, a timestamp-based name will be used.")

    # Paramètre : --deep_analysis
    parser.add_argument("--deep_analysis", type=int, choices=[0, 1], default=0,
                        help="If set to 1, runs extra steps for divergent or novel HPV types. Default: 0.")

    # Paramètre : --plot
    parser.add_argument("--plot", type=str, choices=["pdf", "png"], default="png",
                        help="Specifies the format of the plot (pdf or png). Default is 'png'.")

    return parser.parse_args()


def validate_arguments(args):
    """
    Apply domain validation rules and I/O controls
    to the raw parameters parsed by argparse.
    """
    # 1. Validation du fichier Sample 
    if not Path(args.sample).is_file():
        print(f"Error: Sample file '{args.sample}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # 2. Arbre décisionnel strict pour la base de données (--db)
    db_input = str(args.db).lower() 

    if re.fullmatch(r"pave(_\w+)?", db_input) and not db_input.endswith((".fasta", ".fa", ".fna")):
        # On définit les DEUX nouvelles bases
        db_classic_path = SCRIPT_DIR / "data" / "pave_classic.fasta"
        db_non_official_path = SCRIPT_DIR / "data" / "pave_non_official.fasta"
        print(f"[INFO] Base PAVE (Mode Dual-Search activé). Version : {db_input}")
        args.db = db_input # Maintien de la version (ex: pave_2026) en variable
    else:
        # Validation sécurisée de l'existence du fichier Custom (Correction du bug args.sample)
        if db_input != "genbank" and not Path(args.db).is_file():
            print(f"Error: Database file '{args.db}' does not exist.", file=sys.stderr)
            sys.exit(1)

    # 3. Assainissement du nom de session
    if args.run_name:
        args.run_name = re.sub(r'[^a-zA-Z0-9_-]', '', args.run_name)

    return args



# --- CONFIGURATION NCBI ---
HPV_SEARCH_TERM = (
    '(("Human papillomavirus"[Organism] OR "Human papillomavirus"[All Fields] OR HPV[All Fields]) '
    'AND 7000:10000[Sequence Length]) AND viruses[filter]'
)

def get_remote_genbank_count():
    """Query NCBI to retrieve the number of HPV sequences without downloading them"""
    try:
        with Entrez.esearch(db="nucleotide", term=HPV_SEARCH_TERM, retmax=0) as handle:
            record = Entrez.read(handle)
            return int(record["Count"])
    except Exception as e:
        print(f"[ERROR] Failed to connect to NCBI. : {e}", file=sys.stderr)
        return -1

def download_genbank_fasta(output_path, total_count):
    """Streaming I/O : Batch downloading and on-the-fly filtering."""
    with Entrez.esearch(db="nucleotide", term=HPV_SEARCH_TERM, retmax=total_count) as handle:
        record = Entrez.read(handle)
        id_list = record["IdList"]

    batch_size = 500
    kept_sequences = 0
    dropped_sequences = 0
    
    with open(output_path, "w") as f:
        for start in range(0, len(id_list), batch_size):
            end = min(start + batch_size, len(id_list))
            batch_ids = id_list[start:end]
            
            try:
                logger.info(f"Downloading and filtering sequences {start + 1} to {end}...")
                with Entrez.efetch(db="nucleotide", id=batch_ids, rettype="fasta", retmode="text") as handle:
                    # Lecture du flux NCBI ligne par ligne
                    keep_current = False
                    for line in handle:
                        if line.startswith(">"):
                            # Évaluation C-level du header via l'expression régulière
                            if EXCLUSION_PATTERN.search(line):
                                keep_current = False
                                dropped_sequences += 1
                            else:
                                keep_current = True
                                kept_sequences += 1
                                f.write(line)
                        elif keep_current:
                            # Écriture de l'ADN uniquement si le header a été validé
                            f.write(line)
            except Exception as e:
                logger.error(f"Error processing batch {start}-{end}: {e}")
            
            time.sleep(1) # Respect strict des limites NCBI API
            
    logger.info(f"GenBank filtering complete: {kept_sequences} retained, {dropped_sequences} discarded.")
    return kept_sequences

def setup_genbank_db(args):
    """
    Timestamped local GenBank archive manager.
    Enforces offline execution by default and handles initial cold start.
    """
    if args.db != "genbank":
        return args

    data_dir = SCRIPT_DIR / "data"
    data_dir.mkdir(exist_ok=True)
    
    # Court-circuit : Si l'utilisateur pointe explicitement vers un fichier FASTA spécifique via --db
    default_pave_path = str(data_dir / "DB_full.fasta")
    if str(args.db) != default_pave_path and args.db != "genbank":
        if Path(args.db).is_file():
            logger.info(f"Explicit Archive Mode: NCBI synchronization disabled. Using -> {args.db}")
            return args

    # 1. Balayage du répertoire data/ pour répertorier les archives locales existantes (genbank_YYYY-MM-DD.fasta)
    existing_archives = sorted(list(data_dir.glob("genbank_*.fasta")), reverse=True)
    latest_archive = existing_archives[0] if existing_archives else None

    # 2. Évaluation conditionnelle : Déclenchement réseau
    if args.update_genbank or not latest_archive:
        if not latest_archive:
            logger.info("No local GenBank archive detected. Initializing NCBI download...")
        else:
            logger.info("--update_genbank flag detected. Querying the NCBI Entrez API...")
        
        count_remote = get_remote_genbank_count()
        if count_remote == -1:
            logger.critical("Communication failure with the NCBI Entrez API. Process aborted.")
            sys.exit(1)
            
        today_str = datetime.datetime.now().strftime("%Y-%m-%d")
        new_archive_path = data_dir / f"genbank_{today_str}.fasta"
        
        # Court-circuit 1 : Prévention des téléchargements redondants le même jour calendaire
        if new_archive_path.is_file():
            logger.info(f"A timestamped archive already exists for today's date ({new_archive_path.name}). Skipping download.")
            args.db = str(new_archive_path)
            return args

        # --- NOUVEAU COURT-CIRCUIT 2 : Vérification du delta de séquences ---
        if latest_archive and latest_archive.is_file():
            logger.info(f"Counting sequences in the latest local archive ({latest_archive.name})...")
            # Comptage C-level ultra-rapide des lignes d'en-tête FASTA
            with open(latest_archive, 'rb') as f:
                count_local_archive = sum(1 for line in f if line.startswith(b'>'))
            
            # Si le serveur NCBI n'a pas plus de séquences que notre archive existante, on annule l'I/O réseau
            if count_remote <= count_local_archive:
                logger.info(f"NCBI Entrez ({count_remote}) contains no new sequences compared to the local archive ({count_local_archive}). Download skipped.")
                args.db = str(latest_archive)
                return args
            else:
                logger.info(f"Update required: {count_remote - count_local_archive} new sequences detected on NCBI.")
        # ---------------------------------------------------------------------

        nb_downloaded = download_genbank_fasta(new_archive_path, count_remote)
        logger.info(f"New GenBank archive generated: {new_archive_path.name} ({nb_downloaded} valid sequences).")
        args.db = str(new_archive_path)

        return args
    

def run_post_processing(args, analyze_dir, final_dir, vsearch_classic, vsearch_non_official, out_csv):
    """
    Post-alignment routing architecture.
    Executes the appropriate VSEARCH subprocesses based on the target database.
    """
    genbank_matches_fasta = analyze_dir / "genbank_matches.fasta"
    besthit_genbank_match_tsv = analyze_dir / "besthit_genbank_match.tsv"   

    try:
        # BRANCHE 1 : PAVE 
        if args.is_pave:
            print("[INFO] Primary analysis based on PaVE database")
            excel_file = final_dir / f"{args.run_name}_results.xlsx"
            process_vsearch_clusters(
                vsearch_classic_file=vsearch_classic, 
                vsearch_non_official_file=vsearch_non_official,
                out_csv=out_csv, 
                analyze_dir=analyze_dir, 
                db_fasta=args.db_classic, 
                sample_fasta=args.sample,
                excel_out=excel_file,
                is_genbank=False
            )

        # BRANCHE 2 : GENBANK OU CUSTOM DB
        else:
            is_gb_flag = args.db_classic.name.startswith("genbank_")
            db_type = "GenBank" if is_gb_flag else "Custom"
            print(f"[INFO] Main Analysis Based on {db_type} Database (Cross-aligned with PaVE)")
            
            # Étape A: Parsing initial (L'astuce is_genbank=True s'applique aux deux !)
            process_vsearch_clusters(
                vsearch_classic_file=vsearch_classic,
                vsearch_non_official_file=None,
                out_csv=out_csv, 
                analyze_dir=analyze_dir, 
                db_fasta=args.db_classic, 
                sample_fasta=args.sample,
                is_genbank=True 
            )

            # Étape B: Extraction
            extract_fasta_sequences(
                tsv_path=out_csv,
                fasta_in_path=args.db_classic,
                fasta_out_path=genbank_matches_fasta
            )

            # Étape C: Contre-alignement sur PaVE
            pave_reference = str(SCRIPT_DIR / "data" / "pave_classic_2026.fasta")
            vsearch_cmd_2 = [
                "vsearch", "--usearch_global", str(genbank_matches_fasta),
                "--db", pave_reference, "--id", "0", "--maxaccepts", "1",
                "--threads", "8", "--blast6out", str(besthit_genbank_match_tsv)
            ]
            subprocess.run(vsearch_cmd_2, check=True)

            # Étape D: Jointure
            excel_file = final_dir / f"{args.run_name}_results.xlsx"
            generate_final_custom_report(
                tsv_one_path=out_csv,
                tsv_two_path=besthit_genbank_match_tsv,
                excel_out=excel_file,
                format_mode="compact"
            )

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] VSEARCH subprocess execution failed. Exit code: {e.returncode}", file=sys.stderr)
        sys.exit(1)

def run_deep_analysis(args, analyze_dir, final_dir):
    """
    Execute in-depth L1 gene analysis if divergent types were detected during primary analysis.
    """
    # 1. Évaluation de la première porte logique (L'option utilisateur)
    if args.deep_analysis != 1 or args.db.lower() != "pave":
        return

    # Utilisation stricte de pathlib pour la résolution des chemins
    probable_newtypes_txt = analyze_dir / "probable_newtypes.txt"
    scripts_dir = SCRIPT_DIR / "scripts"

    # 2. Évaluation de la seconde porte logique (Existence du fichier déclencheur)
    if probable_newtypes_txt.is_file():
        print("[INFO] Deep analysis of divergent types...")

        # Préparation des chemins via l'opérateur de surcharge /
        prob_newtypes_fasta = analyze_dir / f"{args.run_name}_prob_newtypes.fasta"
        puma_results_dir = analyze_dir / f"{args.run_name}_puma_results"
        samples_l1_fasta = analyze_dir / f"{args.run_name}_samples_L1.fasta"
        db_full_l1 = SCRIPT_DIR / "data" / "DB_full_L1.fasta"
        sample_l1_besthit = analyze_dir / f"{args.run_name}_sample_L1.besthit.tsv"
        result_deep_xlsx = final_dir / f"{args.run_name}_result_deep_analysis.xlsx"

        try:
            # Étape 1 : Extraction des séquences entières filtrées (Remplace catch_new.pl)
            prepare_deep_analysis_fasta(
                fasta_in_path=args.sample,
                id_list_path=probable_newtypes_txt,
                fasta_out_path=prob_newtypes_fasta
            )

            # Étape 2 : Création du dossier et exécution de PUMA
            # (PUMA est un outil externe écrit en Python, le subprocess est pertinent ici)
            puma_results_dir.mkdir(parents=True, exist_ok=True)
            script_puma = scripts_dir / "puma" / "scripts" / "run_puma.py"
            puma_data_dir = scripts_dir / "puma" / "data_dir"
            subprocess.run(["python3", str(script_puma), "-i", str(prob_newtypes_fasta), 
                            "-d", str(puma_data_dir), "-o", str(puma_results_dir)], check=True)
            # Étape 3 : Isolement du gène L1 (Remplace cut_L1.pl)
            extract_l1_sequences(
                puma_results_dir=puma_results_dir,
                output_fasta=samples_l1_fasta
            )

            # Étape 4 : Ré-alignement ciblé sur L1 (Outil externe, reste en subprocess)
            vsearch_cmd_3 = [
                "vsearch",
                "--usearch_global", str(samples_l1_fasta),
                "--db", str(db_full_l1),
                "--id", "0",
                "--maxaccepts", "1",
                "--threads", "8",
                "--blast6out", str(sample_l1_besthit)
            ]
            subprocess.run(vsearch_cmd_3, check=True)

            # Étape 5 : Calcul des distances et génération du TSV
            finalize_deep_analysis(
                vsearch_l1_file=sample_l1_besthit,
                output_excel=result_deep_xlsx
            )
        except subprocess.CalledProcessError:
            pass

# ==========================================
# ZONE D'EXÉCUTION PRINCIPALE (MAIN)
# ==========================================
if __name__ == "__main__":
    
    # 1. Initialisation et validation
    args_bruts = parse_arguments()
    args_valides = validate_arguments(args_bruts)

    # 2. Synchronisation GenBank (Protégée)
    args_finaux = args_valides
    if args_finaux and "genbank" in str(args_finaux.db).lower():
        args_finaux = setup_genbank_db(args_valides)
        
    # Sécurité ultime au cas où la fonction renverrait quand même du vide
    if args_finaux is None:
        args_finaux = args_valides

    # 3. Nom de session par défaut
    if not args_finaux.run_name:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") 
        args_finaux.run_name = f"run_{timestamp}"

    # 4. Création des répertoires de travail
    run_dir = SCRIPT_DIR / "results" / args_finaux.run_name
    
    analyze_dir = run_dir / "intermediates" 
    final_dir = run_dir / "reports"
    
    analyze_dir.mkdir(parents=True, exist_ok=True)
    final_dir.mkdir(parents=True, exist_ok=True)


    # 5. Définition des fichiers intermédiaires
    out_csv = analyze_dir / f"results_{args_finaux.run_name}.tsv"
    vsearch_out = analyze_dir / f"{args_finaux.run_name}_sample.tsv"

    sanitized_sample = analyze_dir / f"{args_finaux.run_name}_sanitized_sample.fasta"
    try:
        logger.info("Sequence sanitization and removal of invalid nucleotides.")
        sanitize_sample_fasta(args_finaux.sample, sanitized_sample)
    except DataValidationError as e:
        logger.critical(str(e))
        sys.exit(1)

    # 6. Alignement VSEARCH 
    logger.info(f"Starting VSEARCH sequence alignment for run: {args_finaux.run_name}")
    
    # --- DÉFINITION LOCALE DES CHEMINS (AVEC VERSIONING PAVE ET GENBANK) ---
    db_arg_lower = str(args_finaux.db).lower()
    is_pave_mode = bool(re.fullmatch(r"pave(_\w+)?", db_arg_lower)) and not db_arg_lower.endswith((".fasta", ".fa", ".fna"))
    
    
    if is_pave_mode:
        data_dir = SCRIPT_DIR / "data"
        
        # Extraction du suffixe de version saisi par l'utilisateur
        version_suffix = db_arg_lower.replace("pave", "") 
        
        # Assignation de la version par défaut si l'argument strict "--db pave" est utilisé
        if not version_suffix:
            version_suffix = "_2026"
            
        # Construction des chemins absolus
        db_classic_path = data_dir / f"pave_classic{version_suffix}.fasta"
        db_non_official_path = data_dir / f"pave_non_official{version_suffix}.fasta"
        
        # Validation de l'existence des fichiers physiques
        if not db_classic_path.exists() or not db_non_official_path.exists():
            logger.critical(f"Error: PaVE reference files not found ({db_classic_path.name} or {db_non_official_path.name})")
            sys.exit(1)
    else:
        # Résolution dynamique pour GenBank ou Base Custom
        if str(args_finaux.db).lower() == "genbank":
            data_dir = SCRIPT_DIR / "data"
            # Cherche la dernière archive téléchargée (ex: genbank_2026-10-10.fasta)
            archives = sorted(list(data_dir.glob("genbank_*.fasta")), reverse=True)
            db_classic_path = archives[0] if archives else data_dir / "genbank_database.fasta"
        else:
            db_classic_path = Path(args_finaux.db)
            
        db_non_official_path = None
    # -------------------------------------
    try:
        if is_pave_mode:
            # === DUAL-SEARCH PAVE ===
            vsearch_classic = analyze_dir / f"{args_finaux.run_name}_classic.tsv"
            vsearch_non_official = analyze_dir / f"{args_finaux.run_name}_msk.tsv"
            
            vsearch_classic_cmd = [
                "vsearch", "--usearch_global", str(sanitized_sample), "--db", str(db_classic_path),
                "--id", "0", "--maxaccepts", "1", "--threads", "8", "--blast6out", str(vsearch_classic)
            ]
            
            vsearch_non_official_cmd = [
                "vsearch", "--usearch_global", str(sanitized_sample), "--db", str(db_non_official_path),
                "--id", "0", "--maxaccepts", "1", "--threads", "8", "--blast6out", str(vsearch_non_official)
            ]
            
            subprocess.run(vsearch_classic_cmd, check=True)
            subprocess.run(vsearch_non_official_cmd, check=True)
            logger.info("[INFO] PaVE alignment completed successfully.")
        else:
            # === SINGLE-SEARCH GENBANK/CUSTOM ===
            vsearch_classic = analyze_dir / f"{args_finaux.run_name}_sample.tsv"
            vsearch_non_official = None
            
            vsearch_cmd = [
                "vsearch", "--usearch_global", str(sanitized_sample), "--db", str(db_classic_path),
                "--id", "0", "--maxaccepts", "1", "--threads", "8", "--blast6out", str(vsearch_classic)
            ]
            subprocess.run(vsearch_cmd, check=True)
            logger.info("Alignement unique terminé avec succès.")
            
    except FileNotFoundError:
        logger.critical("f[ERROR] 'vsearch' executable not found or not present in the system PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logger.critical(f"[ERROR] VSEARCH execution failed (exit code: {e.returncode}).")
        sys.exit(1)

    # Redirection de l'argument vers le fichier nettoyé pour la suite du pipeline
    args_finaux.sample = str(sanitized_sample)

    # --- ATTACHEMENT DES VARIABLES POUR LE POST-PROCESSING ---
    args_finaux.is_pave = is_pave_mode
    args_finaux.db_classic = db_classic_path
    args_finaux.db_non_official = db_non_official_path
    # ---------------------------------------------------------

    # 7. Post-Processing et Rapports
    run_post_processing(args_finaux, analyze_dir, final_dir, vsearch_classic, vsearch_non_official, out_csv)

    # 8. Analyse approfondie (PUMA / L1)
    run_deep_analysis(args_finaux, analyze_dir, final_dir)

# ======= FINAL PLOTS =======
    if args_finaux.is_pave:
        excel_results = final_dir / f"{args_finaux.run_name}_results.xlsx"
        if excel_results.exists():
            print(f"[INFO] Generating distribution plots ({args_finaux.plot})...")
            generate_plots(
                excel_file=excel_results,
                session_name=args_finaux.run_name,
                output_dir=final_dir,
                fmt=args_finaux.plot,
                mode="compact"
            )
            print(f"[INFO] Plots saved in: {final_dir}")
        else:
            print("[WARNING] PaVE Results Excel file not found. Skipping plot generation.")
    else:
        # Message propre et informatif pour les bases Custom et GenBank
        print("[INFO] Custom/GenBank run detected. Plot generation is exclusive to PaVE (Skipping).")

    print("\n[SUCCÈS] Pipeline execution completed.")