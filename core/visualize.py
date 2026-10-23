import matplotlib       # 1. On importe d'abord le module principal
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from pathlib import Path
import logging

def setup_style():
    """Configuration stricte des paramètres graphiques Matplotlib/Seaborn."""
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'axes.edgecolor': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'font.size': 12
    })

def save_plot(filename, output_dir, fmt, width, height, use_tight_layout=True):
    """Encapsulation de la sauvegarde avec application de la résolution."""
    path = Path(output_dir) / f"{filename}.{fmt}"
    
    # Appel conditionnel pour éviter les conflits de redimensionnement sur les FacetGrids
    if use_tight_layout:
        plt.tight_layout()
        
    if fmt == 'png':
        plt.savefig(path, dpi=300, bbox_inches='tight')
    else:
        plt.savefig(path, bbox_inches='tight')
    plt.close()

def generate_plots(excel_file, session_name, output_dir, fmt="png", mode="compact"):
    """
    Traitement vectoriel et génération des distributions taxonomiques.
    """
    setup_style()
    output_dir = Path(output_dir)
    
    try:
        df = pd.read_excel(excel_file)
    except Exception as e:
        print(f"[ERROR] Échec de la lecture du fichier de résultats : {e}")
        return

    # 1. Extraction et formatage des colonnes cibles
    if 'PAVE match' in df.columns and mode == "compact":
        # Séparation vectorisée avec prévention des dépassements d'index
        temp_split = df['PAVE match'].astype(str).str.split('-', expand=True)
        data = pd.DataFrame({
            'Type': temp_split[0] if temp_split.shape[1] > 0 else pd.NA,
            'Lineage': temp_split[1] if temp_split.shape[1] > 1 else pd.NA,
            'Sub-lineage': temp_split[2] if temp_split.shape[1] > 2 else pd.NA
        })
    else:
        # Correspondance exacte des colonnes attendues
        expected_cols = ['Type', 'Lineage', 'Sub-lineage']
        available_cols = [c for c in expected_cols if c in df.columns]
        data = df[available_cols].copy()
        for col in expected_cols:
            if col not in data.columns:
                data[col] = pd.NA

    # 2. Nettoyage et normalisation
    for col in data.columns:
        data[col] = data[col].replace(['ND', 'NA', 'nan', '', None], pd.NA)
        data[col] = data[col].str.strip()

    data = data.dropna(subset=['Type']) 
    
    # Remplacement de .apply(lambda) par la vectorisation C-level stricte
    data['Type'] = data['Type'].astype(str).str.replace(r'[^A-Za-z0-9_-]', '', regex=True)
    data = data[data['Type'] != '']

    # --- NOUVELLE CONDITION STRICTE ---
    # La condition d'arrêt vérifie maintenant si la taille du jeu de données est <= 1
    if len(data) < 2:
        print(f"[WARNING] Insufficient data ({len(data)} valid strains detected). Distribution plot skipped.")
        return  

    # --- PLOT 1: Distribution des Types ---
    top_20 = data['Type'].value_counts().head(20)
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x=top_20.index, y=top_20.values, color='steelblue')
    plt.xticks(rotation=45, ha='right')
    
    title_1 = f"Top {len(top_20)} HPV Types" if len(top_20) < 20 else "Top 20 Types"
    plt.title(title_1)
    plt.ylabel("Number of Sequences")
    
    # Décalage vertical dynamique basé sur la valeur
    for i, v in enumerate(top_20.values):
        ax.text(i, v + (v * 0.02), str(int(v)), ha='center', size=10)
    save_plot(f"{session_name}_types", output_dir, fmt, 10, 6)


    # --- PLOT 2: Distribution des Lignées par Type ---
    top_types = top_20.index
    lineage_data = data[(data['Type'].isin(top_types)) & (data['Lineage'].notna())]
    
    if not lineage_data.empty:
        lineage_counts = lineage_data.groupby(['Type', 'Lineage']).size().reset_index(name='count')
        
        # NOUVEAU : On filtre l'axe X pour retirer les types sans lignée (les anciens "ND")
        # tout en conservant l'ordre décroissant initial
        valid_types = [t for t in top_types if t in lineage_counts['Type'].values]
        
        plt.figure(figsize=(12, 7))
        # Utilisation de valid_types pour l'argument order
        sns.barplot(data=lineage_counts, x='Type', y='count', hue='Lineage', palette='colorblind', order=valid_types)
        plt.xticks(rotation=45, ha='right')
        
        # Optionnel : Mettre à jour le titre pour refléter le nombre exact affiché
        plt.title(f"Lineages per HPV Type (Top {len(valid_types)})")
        plt.ylabel("Number of Sequences")
        save_plot(f"{session_name}_lineages_per_type", output_dir, fmt, 12, 7)

    # --- PLOT 3: Distribution des Sous-Lignées (Top 10 avec Faceting) ---
    top_10_types = data['Type'].value_counts().head(10).index
    sub_data = data[(data['Type'].isin(top_10_types)) & (data['Lineage'].notna()) & (data['Sub-lineage'].notna())].copy()
    
    if not sub_data.empty:
        # Extraction du composant numérique. expand=False est obligatoire pour éviter un DataFrame 2D.
        sub_data['SublineageGroup'] = sub_data['Sub-lineage'].str.extract(r'(\d+)', expand=False)
        sub_data = sub_data.dropna(subset=['SublineageGroup'])

        if not sub_data.empty:
            num_cols = min(3, len(sub_data['Type'].unique()))
            g = sns.catplot(
                data=sub_data, x='Lineage', kind='count', hue='SublineageGroup',
                col='Type', col_wrap=num_cols, palette='Dark2', height=4, aspect=1.2,
                sharex=False, sharey=False
            )
            g.set_axis_labels("Lineage", "Number of Sequences")
            g.fig.suptitle(f"Top {len(top_10_types)} Types - Sublineages per Lineage", y=1.02)
            
            sns.move_legend(g, "center right", bbox_to_anchor=(1.15, 0.5))
            
            g.fig.subplots_adjust(right=0.85)

            save_plot(f"{session_name}_sublineages_top10_types", output_dir, fmt, 14, 10, use_tight_layout=False)