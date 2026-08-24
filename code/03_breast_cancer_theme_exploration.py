import pandas as pd
from pathlib import Path
TCGA_DIR = Path.cwd() / 'data_raw' / 'TCGA_BRCA'
tcga_file = TCGA_DIR / 'tcga_brca_mutation_cnv_rna_subtypes.csv'
tcga = pd.read_csv(tcga_file)
print('TCGA shape:', tcga.shape)
print('\nFirst 20 columns:')
print(tcga.columns[:20].tolist())
print('\nLast 20 columns:')
print(tcga.columns[-20:].tolist())
print('\nColumn suffix counts:')
suffixes = ['_mutation', '_cnv', '_rna']
for suffix in suffixes:
    cols = [c for c in tcga.columns if c.endswith(suffix)]
    print(suffix, ':', len(cols))
from pathlib import Path
import pandas as pd
TCGA_DIR = Path.cwd() / 'data_raw' / 'TCGA_BRCA'
for file in TCGA_DIR.glob('*.csv'):
    print('\n' + '=' * 80)
    print('FILE:', file.name)
    print('=' * 80)
    df = pd.read_csv(file, nrows=3)
    print('Number of columns:', len(df.columns))
    print('First 15 columns:')
    print(df.columns[:15].tolist())
    print('\nLast 15 columns:')
    print(df.columns[-15:].tolist())
import pandas as pd
from pathlib import Path
TCGA_DIR = Path.cwd() / 'data_raw' / 'TCGA_BRCA'
tcga_file = TCGA_DIR / 'tcga_brca_as_validation.csv'
tcga_cols = pd.read_csv(tcga_file, nrows=0).columns.tolist()
print('Total columns:', len(tcga_cols))
feature_cols = [c for c in tcga_cols if c != 'Cell_line']
rna_cols = [c for c in feature_cols if c.endswith('_rna')]
cnv_cols = [c for c in feature_cols if c.endswith('_cnv')]
mutation_cols = [c for c in feature_cols if c.endswith('_mutation')]
print('RNA features:', len(rna_cols))
print('CNV features:', len(cnv_cols))
print('Mutation features:', len(mutation_cols))

def strip_suffix(cols, suffix):
    return {c[:-len(suffix)] for c in cols if c.endswith(suffix)}
rna_genes = strip_suffix(rna_cols, '_rna')
cnv_genes = strip_suffix(cnv_cols, '_cnv')
mutation_genes = strip_suffix(mutation_cols, '_mutation')
tcga_brca_genes = rna_genes | cnv_genes | mutation_genes
print('Unique RNA genes:', len(rna_genes))
print('Unique CNV genes:', len(cnv_genes))
print('Unique mutation genes:', len(mutation_genes))
print('Total unique TCGA-BRCA genes:', len(tcga_brca_genes))
print('\nFirst 30 genes:')
print(sorted(tcga_brca_genes)[:30])
brca_coverage_rows = []
for db, genes in gene_sets.items():
    overlap = genes & tcga_brca_genes
    resource_coverage = len(overlap) / len(genes)
    brca_coverage = len(overlap) / len(tcga_brca_genes)
    brca_coverage_rows.append({'Database': db, 'Resource_Genes': len(genes), 'TCGA_BRCA_Genes': len(tcga_brca_genes), 'Shared_with_TCGA_BRCA': len(overlap), 'Fraction_of_Resource_in_TCGA_BRCA': resource_coverage, 'Fraction_of_TCGA_BRCA_covered': brca_coverage})
brca_coverage_df = pd.DataFrame(brca_coverage_rows)
brca_coverage_df['Fraction_of_Resource_in_TCGA_BRCA'] = brca_coverage_df['Fraction_of_Resource_in_TCGA_BRCA'].round(4)
brca_coverage_df['Fraction_of_TCGA_BRCA_covered'] = brca_coverage_df['Fraction_of_TCGA_BRCA_covered'].round(4)
print(brca_coverage_df)
brca_coverage_df.to_csv(Path.cwd() / 'results' / 'TCGA_BRCA_gene_coverage_by_resource.csv', index=False)
import numpy as np
import matplotlib.pyplot as plt
plot_df = brca_coverage_df.copy()
order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB_Hallmark']
plot_df = plot_df.set_index('Database').loc[order].reset_index()
x = np.arange(len(order))
width = 0.35
plt.figure(figsize=(10, 6))
plt.bar(x - width / 2, plot_df['Fraction_of_Resource_in_TCGA_BRCA'] * 100, width, label='Resource genes represented in TCGA-BRCA')
plt.bar(x + width / 2, plot_df['Fraction_of_TCGA_BRCA_covered'] * 100, width, label='TCGA-BRCA genes covered by resource')
plt.xticks(x, order, rotation=20)
plt.ylabel('Coverage (%)')
plt.xlabel('Pathway resource')
plt.title('Gene coverage of pathway resources in TCGA-BRCA')
plt.legend()
plt.tight_layout()
plt.savefig(Path.cwd() / 'figures' / 'TCGA_BRCA_gene_coverage_comparison.png', dpi=300, bbox_inches='tight')
plt.show()
unique_pathways = master_df[['Database', 'Pathway_ID', 'Pathway_Name']].drop_duplicates().merge(pathway_sizes[['Database', 'Pathway_ID', 'Gene_Count']], on=['Database', 'Pathway_ID'], how='left')
print('Unique pathways table shape:', unique_pathways.shape)
print(unique_pathways.head())
breast_cancer_themes = {'p53': ['p53'], 'PI3K_AKT': ['pi3k', 'akt', 'mtor'], 'ERBB_HER2': ['erbb', 'her2', 'egfr'], 'ESTROGEN': ['estrogen', 'oestrogen', 'esr1'], 'APOPTOSIS': ['apoptosis', 'apoptotic'], 'DNA_REPAIR': ['dna repair', 'homologous recombination', 'mismatch repair', 'base excision repair', 'nucleotide excision repair', 'fanconi'], 'CELL_CYCLE': ['cell cycle', 'g2m', 'e2f', 'checkpoint'], 'EMT': ['epithelial mesenchymal transition', 'emt']}
print('Themes:', list(breast_cancer_themes.keys()))
theme_rows = []
for theme, keywords in breast_cancer_themes.items():
    mask = unique_pathways['Pathway_Name'].str.lower().fillna('').apply(lambda x: any((k.lower() in x for k in keywords)))
    matched = unique_pathways.loc[mask].copy()
    matched['Theme'] = theme
    theme_rows.append(matched)
theme_pathways_df = pd.concat(theme_rows, ignore_index=True)
theme_pathways_df = theme_pathways_df[['Theme', 'Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Count']].drop_duplicates()
print('Matched pathway rows:', len(theme_pathways_df))
print(theme_pathways_df.head(20))
theme_summary = theme_pathways_df.groupby(['Theme', 'Database']).size().reset_index(name='Matched_Pathway_Count')
print(theme_summary)
theme_pathways_df.to_csv(Path.cwd() / 'results' / 'breast_cancer_theme_matched_pathways.csv', index=False)
theme_summary.to_csv(Path.cwd() / 'results' / 'breast_cancer_theme_pathway_summary.csv', index=False)
print('Saved theme-matched pathway files.')
selected_themes = ['p53', 'PI3K_AKT', 'ERBB_HER2', 'ESTROGEN', 'DNA_REPAIR']
theme_gene_sets = {}
for theme in selected_themes:
    theme_gene_sets[theme] = {}
    matched = theme_pathways_df[theme_pathways_df['Theme'] == theme]
    for db in matched['Database'].unique():
        pathway_ids = matched.loc[matched['Database'] == db, 'Pathway_ID'].unique()
        genes = set(master_df.loc[(master_df['Database'] == db) & master_df['Pathway_ID'].isin(pathway_ids), 'Gene_Symbol'].dropna())
        theme_gene_sets[theme][db] = genes
for theme, dbsets in theme_gene_sets.items():
    print('\n', '=' * 60)
    print(theme)
    print('=' * 60)
    for db, genes in dbsets.items():
        print(db, ':', len(genes), 'genes')
from itertools import combinations
theme_overlap_rows = []
for theme, dbsets in theme_gene_sets.items():
    databases = list(dbsets.keys())
    for db1, db2 in combinations(databases, 2):
        g1 = dbsets[db1]
        g2 = dbsets[db2]
        shared = g1 & g2
        union = g1 | g2
        jaccard = len(shared) / len(union) if len(union) > 0 else 0
        theme_overlap_rows.append({'Theme': theme, 'Database_1': db1, 'Database_2': db2, 'Genes_DB1': len(g1), 'Genes_DB2': len(g2), 'Shared_Genes': len(shared), 'Jaccard': jaccard})
theme_overlap_df = pd.DataFrame(theme_overlap_rows)
theme_overlap_df['Jaccard'] = theme_overlap_df['Jaccard'].round(4)
print(theme_overlap_df)
theme_overlap_df.to_csv(Path.cwd() / 'results' / 'breast_cancer_theme_pairwise_gene_overlap.csv', index=False)
core_theme_rows = []
for theme, dbsets in theme_gene_sets.items():
    if len(dbsets) < 2:
        continue
    common_genes = set.intersection(*dbsets.values())
    print(theme, ':', len(common_genes), 'genes common across all available resources')
    for gene in sorted(common_genes):
        core_theme_rows.append({'Theme': theme, 'Gene_Symbol': gene, 'Resource_Count': len(dbsets)})
core_theme_df = pd.DataFrame(core_theme_rows)
core_theme_df.to_csv(Path.cwd() / 'results' / 'breast_cancer_theme_core_genes.csv', index=False)
theme_specific_rows = []
for theme, dbsets in theme_gene_sets.items():
    for db, genes in dbsets.items():
        other_sets = [gset for other_db, gset in dbsets.items() if other_db != db]
        if not other_sets:
            continue
        other_union = set().union(*other_sets)
        specific = genes - other_union
        for gene in sorted(specific):
            theme_specific_rows.append({'Theme': theme, 'Database': db, 'Gene_Symbol': gene})
theme_specific_df = pd.DataFrame(theme_specific_rows)
print(theme_specific_df.groupby(['Theme', 'Database']).size().reset_index(name='Specific_Gene_Count'))
theme_specific_df.to_csv(Path.cwd() / 'results' / 'breast_cancer_theme_resource_specific_genes.csv', index=False)
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
theme = 'p53'
dbsets = theme_gene_sets[theme]
theme_order = [db for db in ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB_Hallmark'] if db in dbsets]
matrix = pd.DataFrame(np.eye(len(theme_order)), index=theme_order, columns=theme_order)
for db1, db2 in combinations(theme_order, 2):
    g1 = dbsets[db1]
    g2 = dbsets[db2]
    j = len(g1 & g2) / len(g1 | g2)
    matrix.loc[db1, db2] = j
    matrix.loc[db2, db1] = j
plt.figure(figsize=(7, 6))
im = plt.imshow(matrix.values, aspect='auto')
plt.colorbar(im, label='Jaccard similarity')
plt.xticks(range(len(theme_order)), theme_order, rotation=30, ha='right')
plt.yticks(range(len(theme_order)), theme_order)
for i in range(len(theme_order)):
    for j in range(len(theme_order)):
        plt.text(j, i, f'{matrix.iloc[i, j]:.2f}', ha='center', va='center')
plt.title('p53-related pathway gene similarity')
plt.tight_layout()
plt.savefig(Path.cwd() / 'figures' / 'p53_pathway_gene_similarity_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()
