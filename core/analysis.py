import pandas as pd
import re
from pathlib import Path
import csv
import logging
import sys

logger = logging.getLogger("EvolutionaryPipeline")

# --- DÉCLARATION DE L'EXCEPTION ---
class DataValidationError(Exception):
    """Raised when input sequence data is corrupted (e.g., presence of ambiguous 'N' nucleotides)."""
    pass

# REGEX format PAVE (HPV16-A1_K02718) HPV[Type]-[Lignée][Sous-lignée]_[Accession]
REGEX_PAVE = re.compile(r'HPV-?([a-zA-Z0-9]+)(?:-([A-Z])(\d*))?', re.IGNORECASE)

HR_TYPES = {"HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV39", "HPV45", 
            "HPV51", "HPV52", "HPV56", "HPV58", "HPV59", "HPV66", "HPV68"}

def get_fasta_lengths(fasta_path):
    """
    Extracts sequence lengths from a FASTA file.
    """
    lengths = {}
    with open(fasta_path, 'r') as f:
        header = None
        seq_len = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    lengths[header] = seq_len
                header = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
        if header:
            lengths[header] = seq_len
            
    return lengths

def process_vsearch_clusters(vsearch_classic_file, vsearch_non_official_file, out_csv, analyze_dir, db_fasta, sample_fasta, excel_out=None, is_genbank=False):
    """
    Parse VSEARCH tabular output and classify nucleotide sequences.
    Dispatches processing logic (PaVE vs. GenBank) according to the is_genbank parameter.
    """
    sample_lengths = get_fasta_lengths(sample_fasta)
    db_lengths = get_fasta_lengths(db_fasta)

    # ==========================================
    # PARSING ET FUSION 
    # ==========================================
    try:
        # 1. Chargement du fichier classique 
        df_classic = pd.read_csv(vsearch_classic_file, sep='\t', header=None, usecols=[0, 1, 2], names=['query', 'target', 'identity'])
        df = df_classic.drop_duplicates(subset=['query'], keep='first').copy()
        
        # 2. Ajout conditionnel des non-officiels (Mode PAVE)
        if not is_genbank and vsearch_non_official_file:
            try:
                df_msk = pd.read_csv(vsearch_non_official_file, sep='\t', header=None, usecols=[0, 1, 2], names=['query', 'msk_target', 'msk_identity'])
                best_msk = df_msk.drop_duplicates(subset=['query'], keep='first').copy()
                
                # Fusion (Left Join) des valeurs brutes pour comparaison
                df = pd.merge(df, best_msk[['query', 'msk_target', 'msk_identity']], on='query', how='left')
                
                # Initialisation de la colonne avec la valeur par défaut
                df['unofficial_match'] = "-"
                
                # Masque booléen : vérifie si l'identité du variant est strictement supérieure à celle du classique
                mask = df['msk_identity'] > df['identity']
                
                # Application du formatage uniquement sur les lignes respectant la condition
                df.loc[mask, 'unofficial_match'] = df.loc[mask, 'msk_target'] + " (" + df.loc[mask, 'msk_identity'].astype(str) + "%)"
                
                # Nettoyage des colonnes temporaires issues de la fusion
                df = df.drop(columns=['msk_target', 'msk_identity'])
                
            except pd.errors.EmptyDataError:
                df['unofficial_match'] = "-"
        else:
            df['unofficial_match'] = "NA"

    except pd.errors.EmptyDataError:
        print("[WARNING] Le fichier VSEARCH classique est vide.")
        cols = ["Sequence", "Match", "Type", "Lineage", "Sub-lineage", "Similarity", "Length Difference", "Results"] if is_genbank else ["Sequence", "Type", "Lineage", "Sub-lineage", "Similarity", "Length Difference", "Results", "Meilleur match non-officiel"]
        pd.DataFrame(columns=cols).to_csv(out_csv, sep='\t', index=False)
        return

    results = []
    prob_newtypes = []

    # ==========================================
    # 2. ITÉRATION ET CLASSIFICATION
    # ==========================================
    for row in df.itertuples(index=False):
        seq_id = str(row.query)
        hit = str(row.target) if pd.notna(row.target) else ''
        perc = float(row.identity) if pd.notna(row.identity) else 0.0
        
        # Extraction sécurisée de la colonne alternative
        unoff_match = str(row.unofficial_match) if pd.notna(row.unofficial_match) else "-"

        # Contournement de l'ajout si le pourçentage est 0.0 (Cas d'un outer join sans hit classique)
        if 0.0 < perc < 95.0:
            prob_newtypes.append(seq_id)

        # ROUTAGE LOGIQUE : GENBANK vs PAVE
        if is_genbank:
            if perc >= 99.5:
                type_val, lineage, sub_lineage, res_str = 'Same type', 'Same lineage', 'Same sublineage', 'Perfect or near-perfect match'
            elif perc >= 98.5:
                type_val, lineage, sub_lineage, res_str = 'Same type', 'Different lineage', 'Different sublineage', 'Close match - minor divergence'
            elif perc >= 98.0:
                type_val, lineage, sub_lineage, res_str = 'Same type', 'Different lineage', 'Different sublineage', 'Minor divergent - review recommended'
            elif perc >= 95.0:
                type_val, lineage, sub_lineage, res_str = 'Same type', 'Different lineage', 'Different sublineage', 'Moderately divergent - review recommended'
            elif perc >=90.0:
                type_val, lineage, sub_lineage, res_str = 'Same type', 'Different lineage', 'Different sublineage', 'Highly divergent - review recommended'
            else:
                type_val, lineage, sub_lineage, res_str = 'ND', 'ND', 'ND', f"--deep_analysis 1 recommended - no strong match ({perc}%)"
            
        else:
            type_val, lineage, sub_lineage = 'ND', 'ND', 'ND'
            
            match_pave = REGEX_PAVE.search(hit)

            if match_pave:
                type_val = f"HPV{match_pave.group(1)}"
                if match_pave.group(2):
                    lineage = match_pave.group(2)
                    sub = match_pave.group(3)
                    sub_lineage = f"{lineage}-{sub}" if sub else f"{lineage}-1"
                else:
                    lineage, sub_lineage = 'A', 'A-1'

            type_ref = type_val

            if perc <= 85.0:
                lineage, sub_lineage = 'ND', 'ND'
                res_str = f"No significant match, try --deep_analysis, the closest type is {type_ref}"
            elif perc >= 99.5:
                res_str = 'Same sublineage'
            elif perc >= 99.0:
                sub_lineage = 'ND'
                res_str = 'Probably new sublineage'
            elif perc >= 98.0:
                lineage, sub_lineage = 'ND', 'ND'
                res_str = 'Probably new lineage - Minor divergent'
            elif perc >= 95.0:
                lineage, sub_lineage = 'ND', 'ND'
                res_str = 'Probably new lineage - Moderately divergent'
            elif perc >= 90.0:
                lineage, sub_lineage = 'ND', 'ND'
                res_str = 'Probably new lineage - Highly divergent'
            else:
                lineage, sub_lineage = 'ND', 'ND'
                res_str = f"No significant match, try --deep_analysis, the closest type is {type_ref}"

            if perc < 90.0:
                type_val = 'ND'

        # CALCUL DES LONGUEURS
        seq_len = sample_lengths.get(seq_id)
        if seq_len is None:
            seq_len = next((v for k, v in sample_lengths.items() if seq_id in k), 0)

        hit_len = db_lengths.get(hit)
        if hit_len is None:
            hit_len = next((v for k, v in db_lengths.items() if hit in k), 0)

        diff = abs(seq_len - hit_len)

        # ==========================================
        # 3. CRÉATION DU DICTIONNAIRE FINAL
        # ==========================================
        
        display_similarity = "ND" if lineage == "ND" else f"{perc}%"
        row_dict = {'Sequence': seq_id}
        if is_genbank:
            row_dict['Match'] = hit
            
        row_dict.update({
            'Type': type_val,
            'Lineage': lineage,
            'Sub-lineage': sub_lineage,
            'Similarity': display_similarity,
            'Length Difference': diff,
            'Results': res_str
        })
        
        # Injection conditionnelle de la nouvelle colonne pour PAVE
        if not is_genbank:
            row_dict['Top Non-Official Match'] = unoff_match
            
        results.append(row_dict)

    # ==========================================
    # 4. GÉNÉRATION DES SORTIES
    # ==========================================
    out_df = pd.DataFrame(results)
    out_df.to_csv(out_csv, sep='\t', index=False)

    if prob_newtypes:
        prob_path = Path(analyze_dir) / "probable_newtypes.txt"
        with open(prob_path, "w") as f:
            for pt in prob_newtypes:
                f.write(f"{pt}\n")

    if excel_out and not is_genbank:
        writer = pd.ExcelWriter(excel_out, engine='xlsxwriter')
        out_df.to_excel(writer, sheet_name='Results', index=False)
        workbook  = writer.book
        worksheet = writer.sheets['Results']

        header_format = workbook.add_format({'bg_color': 'yellow', 'bold': True, 'border': 1})
        hr_format     = workbook.add_format({'bg_color': '#FFC7CE'})

        for col_num, value in enumerate(out_df.columns.values):
            worksheet.write(0, col_num, value, header_format)

        # L'indexation par df.itertuples() gère nativement le changement de dimension du DataFrame
        for i, row_data in enumerate(out_df.itertuples(), start=1):
            if row_data.Type in HR_TYPES:
                worksheet.set_row(i, None, hr_format)

        writer.close()
        print(f"[INFO] Rapport Excel généré : {excel_out}")

def generate_final_custom_report(tsv_one_path, tsv_two_path, excel_out, format_mode="compact"):
    """
    Perform a relational join between Custom and PaVE alignment results,
    then generate an Excel export with conditional formatting.
    """
    # 1. Parsing du fichier VSEARCH contre PAVE (tsv_two)
    try:
        df2 = pd.read_csv(tsv_two_path, sep='\t', header=None, usecols=[0, 1, 2], names=['seq_id', 'hit', 'perc'])
    except pd.errors.EmptyDataError:
        df2 = pd.DataFrame(columns=['seq_id', 'hit', 'perc'])

    pave_data = []
    
    # Extraction taxonomique sur les hits PAVE
    for row in df2.itertuples(index=False):
        seq_id = str(row.seq_id)
        hit = str(row.hit) if pd.notna(row.hit) else ''
        perc = float(row.perc) if pd.notna(row.perc) else 0.0

        type_val, lineage, sub_lineage = 'ND', 'ND', 'ND'

        match_pave = REGEX_PAVE.search(hit)

        if match_pave:
            type_val = f"HPV{match_pave.group(1)}"
            
            if match_pave.group(2):
                lineage = match_pave.group(2)
                sub = match_pave.group(3)
                sub_lineage = f"{lineage}{sub}" if sub else f"{lineage}1"
            else:
                # Fallback strict
                lineage, sub_lineage = 'A', 'A1'

        type_ref = type_val

        # Algorithme de similarité spécifique au format PAVE
        if perc <= 85.0:
            ls, sls, res = 'ND', 'ND', f"No significant match ({perc}), try --deep_analysis, closest type: {type_ref}"
        elif perc >= 99.5:
            ls, sls, res = lineage, sub_lineage, 'Same sublineage'
        elif perc >= 98.5:
            ls, sls, res = lineage, 'ND', 'Probably new sublineage'
        elif perc >= 98.0:
            ls, sls, res = lineage, 'ND', 'Probably new lineage'
        else:
            ls, sls, res = 'ND', 'ND', f"No significant match ({perc}), try --deep_analysis, closest type: {type_ref}"

        # Note technique: le script Perl écrase la différence de longueur par 100 - identité
        diff = round(100.0 - perc, 2)

        if format_mode == "compact":
            pave_match = f"{type_val}-{ls}-{sls}"
            pave_data.append({'Match_ID': seq_id, 'PAVE match': pave_match})
        else:
            pave_data.append({
                'Match_ID': seq_id,
                'Sequence blasted in PAVE': seq_id,
                'PAVE_Type': type_val,
                'PAVE_Lineage': ls,
                'PAVE_Sub-lineage': sls,
                'Percentage': perc,
                'PAVE_Length_Difference': diff,
                'PAVE result': res
            })

    df_pave = pd.DataFrame(pave_data)

    # 2. Parsing du premier TSV (Custom/GenBank)
    df1 = pd.read_csv(tsv_one_path, sep='\t')
    match_col = 'Match' if 'Match' in df1.columns else df1.columns[1]

    # 3. Jointure à gauche (Left Join)
    if df_pave.empty:
        df_final = df1.copy()
        if format_mode == "compact":
            df_final['PAVE match'] = "NA"
        else:
            cols_to_add = ['Sequence blasted in PAVE', 'Type', 'Lineage', 'Sub-lineage', 'Percentage', 'Length Difference', 'PAVE result']
            for col in cols_to_add:
                df_final[col] = "NA"
    else:
        df_final = pd.merge(df1, df_pave, left_on=match_col, right_on='Match_ID', how='left')
        df_final.drop(columns=['Match_ID'], inplace=True)
        df_final.fillna("NA", inplace=True)

        if format_mode == "detailed":
            df_final.rename(columns={
                'PAVE_Type': 'Type',
                'PAVE_Lineage': 'Lineage',
                'PAVE_Sub-lineage': 'Sub-lineage',
                'PAVE_Length_Difference': 'Length Difference'
            }, inplace=True)

    # 4. Écriture via XlsxWriter
    writer = pd.ExcelWriter(excel_out, engine='xlsxwriter')
    df_final.to_excel(writer, sheet_name='Results', index=False)
    
    workbook  = writer.book
    worksheet = writer.sheets['Results']

    header_format = workbook.add_format({'bg_color': 'yellow', 'bold': True, 'border': 1})
    hr_format     = workbook.add_format({'bg_color': '#FFC7CE'})

    for col_num, value in enumerate(df_final.columns.values):
        worksheet.write(0, col_num, value, header_format)

    # Détermination de l'index de colonne pour le formatage conditionnel HR
    target_col_idx = len(df_final.columns) - 1 if format_mode == "compact" else 9

    for row_idx, row_data in enumerate(df_final.itertuples(index=False), start=1):
        cell_value = str(row_data[target_col_idx])
        
        # Regex pour intercepter la nomenclature HPV et vérifier dans le set HR_TYPES
        match = re.search(r'(HPV\d+)', cell_value)
        if match and match.group(1) in HR_TYPES:
            worksheet.write(row_idx, target_col_idx, cell_value, hr_format)

    writer.close()
    print(f"[INFO] Excel relational report generated: {excel_out}")


def extract_fasta_sequences(tsv_path, fasta_in_path, fasta_out_path):
    """
    Filtre un fichier FASTA de référence pour n'extraire que les séquences
    dont les identifiants sont listés dans la colonne 'Match' du TSV.
    Remplace extract_genbank_matches.pl.
    """
    # 1. Récupération des identifiants cibles
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        if df.empty or len(df.columns) < 2:
            print("[WARNING] Empty or malformed TSV file. No sequences to extract.")
            open(fasta_out_path, 'w').close() # Crée un fichier vide
            return
            
        # La colonne 1 (index 1) contient les "Matches" GenBank/Custom
        match_col = df.columns[1]
        
        # On utilise un set() pour éliminer les doublons et avoir une recherche en O(1)
        wanted_ids = set(df[match_col].dropna().astype(str))
        
    except FileNotFoundError:
        print("[WARNING] Empty or malformed TSV file. No sequences to extract.")
        return

    # 2. Filtrage du FASTA (lecture en flux)
    extracted_count = 0
    with open(fasta_in_path, 'r') as f_in, open(fasta_out_path, 'w') as f_out:
        write_flag = False
        for line in f_in:
            if line.startswith('>'):
                # Extraction de l'ID (sans le '>' et coupé au premier espace)
                seq_id = line[1:].split()[0]
                write_flag = (seq_id in wanted_ids)
                if write_flag:
                    extracted_count += 1
            
            if write_flag:
                f_out.write(line)
                
    print(f"[INFO] FASTA extraction completed: {extracted_count} reference sequences isolated to {fasta_out_path}")

def prepare_deep_analysis_fasta(fasta_in_path, id_list_path, fasta_out_path):
    """
    Extrait les séquences listées dans probable_newtypes.txt.
    Applique un filtrage sur la longueur minimale (200 pb) et la pureté nucléotidique (ACGT).
    Formate l'en-tête de sortie au format >ID|ID pour la compatibilité PUMA.
    Remplace catch_new.pl.
    """
    # 1. Chargement des identifiants cibles (recherche en temps constant)
    try:
        with open(id_list_path, 'r') as f:
            wanted_ids = {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        print(f"[ERROR] Identifier list file not found: {id_list_path}")
        return

    if not wanted_ids:
        print("[INFO] No identifiers found. Creating an empty file.")
        open(fasta_out_path, 'w').close()
        return

    # Ensemble des caractères autorisés pour l'évaluation de pureté
    valid_bases = set('ACGTacgt')
    extracted_count = 0

    # 2. Parsing FASTA, filtrage et écriture
    with open(fasta_in_path, 'r') as f_in, open(fasta_out_path, 'w') as f_out:
        header = ''
        seq_acc = []
        
        def process_sequence(current_header, current_seq_acc):
            nonlocal extracted_count
            if current_header and current_header in wanted_ids:
                full_seq = "".join(current_seq_acc)
                # Application des critères de qualité
                if len(full_seq) >= 200 and set(full_seq).issubset(valid_bases):
                    
                    # --- AJOUT : Nettoyage strict de l'en-tête pour PUMA ---
                    # Création d'une table de traduction (Espace, Parenthèse gauche, Parenthèse droite -> Underscores)
                    trans_table = str.maketrans(" ()", "___")
                    clean_header = current_header.translate(trans_table)
                    # -------------------------------------------------------

                    # Formatage strict >ID|ID exigé par PUMA avec l'en-tête propre
                    f_out.write(f">{clean_header}|{clean_header}\n{full_seq}\n")
                    extracted_count += 1

        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                process_sequence(header, seq_acc)
                # Extraction de l'ID en isolant le premier bloc après '>'
                header = line[1:].split()[0]
                seq_acc = []
            else:
                seq_acc.append(line)
        
        # Traitement de la dernière séquence du buffer
        process_sequence(header, seq_acc)

    print(f"[INFO] PUMA preparation completed: {extracted_count} sequences extracted to {fasta_out_path}")

def extract_l1_sequences(puma_results_dir, output_fasta):
    """
    Parcourt récursivement le dossier des résultats PUMA pour localiser les fichiers CSV.
    Extrait la séquence du gène L1 (4ème colonne) et génère un fichier FASTA aggrégé.
    Remplace cut_L1.pl.
    """
    puma_dir = Path(puma_results_dir)
    extracted_count = 0
    
    # Mode 'w' : écrase le fichier s'il existe déjà
    with open(output_fasta, 'w') as f_out:
        # rglob('*.csv') trouve tous les fichiers CSV dans le dossier et ses sous-dossiers
        for csv_file in puma_dir.rglob('*.csv'):
            try:
                with open(csv_file, 'r', newline='', encoding='utf-8') as f_in:
                    reader = csv.reader(f_in)
                    for row in reader:
                        # Vérification de l'existence des colonnes et de la correspondance L1
                        if len(row) >= 4 and row[1].strip() == "L1":
                            # csv_file.stem récupère le nom du fichier sans l'extension
                            f_out.write(f">{csv_file.stem}\n{row[3].strip()}\n")
                            extracted_count += 1
                            break  # Stoppe la lecture de ce fichier dès que L1 est trouvé
            except Exception as e:
                print(f"[WARNING] Read error on file {csv_file.name}: {e}")
                
    print(f"[INFO] L1 extraction completed: {extracted_count} sequences isolated to {output_fasta}")

def finalize_deep_analysis(vsearch_l1_file, output_excel):
    """
    Analyse les résultats VSEARCH du gène L1 et génère le rapport Excel final.
    Remplace find_newtype.pl et la conversion Pandas manuelle dans le main.
    """
    try:
        # --- MODIFICATION 1 : Ajout de la colonne 1 ('hit') dans la lecture ---
        # 0 = Query (Séquence), 1 = Target (Hit), 2 = Identity (%)
        df = pd.read_csv(vsearch_l1_file, sep='\t', header=None, usecols=[0, 1, 2],
                         names=['Sequence', 'hit', 'Identity'])
    except pd.errors.EmptyDataError:
        logger.warning("[INFO] No L1 results to analyze.") # (Ou print selon votre standard)
        return

    # --- MODIFICATION 2 : Extraction du "Plus Proche Voisin" ---
    def extract_closest_type(hit_string):
        match = REGEX_PAVE.search(str(hit_string))
        return f"HPV{match.group(1)}" if match else "ND"

    df['Closest Known Type'] = df['hit'].apply(extract_closest_type)
    # -----------------------------------------------------------

    # Application de la règle taxonomique (Seuil 90%)
    df['Puma result'] = df['Identity'].apply(
        lambda x: "is a new HPV type" if x < 90.0 else "is NOT a new HPV type"
    )

    # Renommage de la colonne pour correspondre à l'en-tête souhaité
    df.rename(columns={'Sequence': 'Sequence with probable new HPV type'}, inplace=True)

    # --- MODIFICATION 3 : Intégration dans le tableau final ---
    # On ajoute l'Identity et le Closest Known Type pour que le tableau soit compréhensible
    final_df = df[['Sequence with probable new HPV type', 'Identity', 'Closest Known Type', 'Puma result']]

    # Sauvegarde directe en Excel
    final_df.to_excel(output_excel, index=False)
    logger.info(f"[INFO] Final in-depth analysis report generated: {output_excel}")

def sanitize_sample_fasta(input_fasta_path, sanitized_fasta_path, min_len=7000, max_len=10000):
    """
    I/O Streaming : Convertit les séquences en majuscules, supprime celles
    contenant des caractères hors de {A, C, G, T} et filtre par taille.
    Complexité mémoire : O(1) (Une séquence en RAM à la fois).
    """
    valid_bases = {"A", "C", "G", "T"}
    kept_count = 0
    dropped_chars = []
    dropped_length = []

    with open(input_fasta_path, 'r') as f_in, open(sanitized_fasta_path, 'w') as f_out:
        current_header = ""
        current_seq_acc = []

        def process_record():
            nonlocal kept_count
            if not current_header:
                return
            
            # Reconstruction et majuscule
            full_seq = "".join(current_seq_acc).upper()
            seq_len = len(full_seq)
            
            # 1. Évaluation C-level de l'appartenance à l'alphabet
            if not set(full_seq).issubset(valid_bases) or seq_len == 0:
                invalid_chars = set(full_seq) - valid_bases
                dropped_chars.append((current_header, invalid_chars))
                return
            
            # 2. NOUVEAU : Évaluation des limites de taille
            if not (min_len <= seq_len <= max_len):
                dropped_length.append((current_header, seq_len))
                return
                
            # Si la séquence passe les deux tests, on l'écrit
            f_out.write(f"{current_header}\n{full_seq}\n")
            kept_count += 1

        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith(">"):
                process_record()
                current_header = line
                current_seq_acc = []
            else:
                current_seq_acc.append(line)
        
        process_record() # Traitement du dernier buffer

    # Observabilité (Remplacement de logger par print)
    if dropped_chars:
        print(f"[WARNING] Strict filtering: Rejected {len(dropped_chars)} sequence(s) containing non-conforming characters (e.g., 'N').", file=sys.stderr)
    
    if dropped_length:
        print(f"[WARNING] Length filtering: Rejected {len(dropped_length)} sequence(s) outside thresholds ({min_len}-{max_len} bp).", file=sys.stderr)

    if kept_count == 0:
        raise DataValidationError(
            f"Invalid input FASTA '{input_fasta_path}': All sequences were discarded "
            f"after strict filtering (unrecognized characters or length boundary violations)."
        )

    print(f"[INFO] Sanitization complete: {kept_count} sequence(s) retained in {sanitized_fasta_path}.")
    return kept_count