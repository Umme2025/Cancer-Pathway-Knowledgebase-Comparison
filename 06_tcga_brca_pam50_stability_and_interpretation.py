from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
TCGA_FILE = PROJECT_ROOT / 'data_raw' / 'TCGA_BRCA' / 'tcga_brca_as_validation.csv'
tcga = pd.read_csv(TCGA_FILE, nrows=5)
print('TCGA shape preview:', tcga.shape)
print('\nFirst 15 columns:')
print(tcga.columns[:15].tolist())
print('\nPossible ID / label columns:')
for c in tcga.columns:
    name = str(c).lower()
    if any((x in name for x in ['cell', 'sample', 'patient', 'subtype', 'pam50', 'label', 'cancer_type'])):
        print(c)
print('\nPossible subtype/label files in project:')
for p in PROJECT_ROOT.rglob('*.csv'):
    n = p.name.lower()
    if any((x in n for x in ['subtype', 'pam50', 'label', 'tcga_brca'])):
        print(p)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
print('Possible subtype files:\n')
for p in PROJECT_ROOT.rglob('*.csv'):
    name = p.name.lower()
    if any((x in name for x in ['subtype', 'pam50', 'label', 'brca'])):
        try:
            preview = pd.read_csv(p, nrows=3)
            print('\nFILE:', p)
            print('Columns:', preview.columns.tolist()[:20])
        except Exception:
            pass
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
TCGA_FEATURE_FILE = PROJECT_ROOT / 'data_raw' / 'TCGA_BRCA' / 'tcga_brca_as_validation.csv'
TCGA_LABEL_FILE = PROJECT_ROOT / 'data_raw' / 'TCGA_BRCA' / 'tcga_brca_mutation_cnv_rna_subtypes.csv'
RESULTS_DIR = PROJECT_ROOT / 'results'
features = pd.read_csv(TCGA_FEATURE_FILE)
labels = pd.read_csv(TCGA_LABEL_FILE)
print('Feature shape:', features.shape)
print('Label shape:', labels.shape)
print('\nLabel columns:')
print(labels.columns.tolist())
labels = labels.rename(columns={'Cancer_type': 'PAM50_Subtype'})
labels['Cell_line'] = labels['Cell_line'].astype(str).str.strip()
features['Cell_line'] = features['Cell_line'].astype(str).str.strip()
print('\nDuplicate feature IDs:', features['Cell_line'].duplicated().sum())
print('Duplicate label IDs:', labels['Cell_line'].duplicated().sum())
merged = features.merge(labels[['Cell_line', 'PAM50_Subtype']], on='Cell_line', how='left')
matched = merged['PAM50_Subtype'].notna().sum()
unmatched = merged['PAM50_Subtype'].isna().sum()
print('\nMerged shape:', merged.shape)
print('Matched samples:', matched)
print('Unmatched samples:', unmatched)
print('\nSubtype counts:')
print(merged['PAM50_Subtype'].value_counts(dropna=False))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
merged.to_csv(OUT_FILE, index=False)
print('\nSaved:')
print(OUT_FILE)
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kruskal
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
GENE_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
tcga = pd.read_csv(TCGA_FILE)
genes = pd.read_csv(GENE_FILE)
print('TCGA:', tcga.shape)
print('Family-gene table:', genes.shape)
consensus = genes[genes['Resource_Support'] >= 2].copy()
consensus['Gene_Symbol'] = consensus['Gene_Symbol'].astype(str).str.upper().str.strip()
print('Consensus family-gene rows:', len(consensus))
needed_genes = set(consensus['Gene_Symbol'])
rna_cols = [c for c in tcga.columns if c.lower().endswith('_rna') and c[:-4].upper() in needed_genes]
cnv_cols = [c for c in tcga.columns if c.lower().endswith('_cnv') and c[:-4].upper() in needed_genes]
feature_cols = rna_cols + cnv_cols
print('RNA features used:', len(rna_cols))
print('CNV features used:', len(cnv_cols))
X = tcga[feature_cols].apply(pd.to_numeric, errors='coerce')
means = X.mean(axis=0)
stds = X.std(axis=0, ddof=0).replace(0, np.nan)
Z = (X - means) / stds
score_table = tcga[['Cell_line', 'PAM50_Subtype']].copy()
family_info = []
for family_id, g in consensus.groupby('Family_ID'):
    family_genes = set(g['Gene_Symbol'])
    fam_rna = [c for c in rna_cols if c[:-4].upper() in family_genes]
    fam_cnv = [c for c in cnv_cols if c[:-4].upper() in family_genes]
    rna_score = Z[fam_rna].mean(axis=1) if fam_rna else pd.Series(np.nan, index=Z.index)
    cnv_score = Z[fam_cnv].mean(axis=1) if fam_cnv else pd.Series(np.nan, index=Z.index)
    combined_score = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1)
    score_table[f'{family_id}_RNA'] = rna_score
    score_table[f'{family_id}_CNV'] = cnv_score
    score_table[f'{family_id}_Combined'] = combined_score
    family_info.append({'Family_ID': family_id, 'Consensus_Genes': len(family_genes), 'RNA_Features': len(fam_rna), 'CNV_Features': len(fam_cnv)})
family_info = pd.DataFrame(family_info)
subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
results = []
for family_id in family_info['Family_ID']:
    col = f'{family_id}_Combined'
    groups = [score_table.loc[score_table['PAM50_Subtype'] == subtype, col].dropna().values for subtype in subtypes]
    H, p = kruskal(*groups)
    n = sum((len(x) for x in groups))
    k = len(groups)
    epsilon2 = max(0, (H - k + 1) / (n - k))
    row = {'Family_ID': family_id, 'Kruskal_H': H, 'P_Value': p, 'Epsilon_Squared': epsilon2}
    for subtype in subtypes:
        row[f'Median_{subtype}'] = score_table.loc[score_table['PAM50_Subtype'] == subtype, col].median()
    results.append(row)
association = pd.DataFrame(results)
association = association.sort_values('P_Value').reset_index(drop=True)
m = len(association)
association['FDR'] = association['P_Value'] * m / np.arange(1, m + 1)
association['FDR'] = association['FDR'][::-1].cummin()[::-1].clip(upper=1)
association['Significant_FDR_0.05'] = association['FDR'] < 0.05
association = association.merge(family_info, on='Family_ID', how='left')
print('\n' + '=' * 100)
print('STEP 9B PAM50 FAMILY ASSOCIATION')
print('=' * 100)
print('\nSignificant families (FDR < 0.05):', association['Significant_FDR_0.05'].sum(), '/', len(association))
display(association[['Family_ID', 'Kruskal_H', 'P_Value', 'FDR', 'Epsilon_Squared', 'Significant_FDR_0.05', 'Median_Basal', 'Median_Her2', 'Median_LumA', 'Median_LumB', 'Median_Normal-like', 'Consensus_Genes']].round(5))
SCORE_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
ASSOC_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
INFO_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_family_features_used.csv'
score_table.to_csv(SCORE_FILE, index=False)
association.to_csv(ASSOC_FILE, index=False)
family_info.to_csv(INFO_FILE, index=False)
print('\nSTEP 9B COMPLETE')
print('\nSaved:')
print(SCORE_FILE)
print(ASSOC_FILE)
print(INFO_FILE)
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
ASSOC_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
association = pd.read_csv(ASSOC_FILE)
top10 = association.sort_values('Epsilon_Squared', ascending=False).head(10).sort_values('Epsilon_Squared')
plt.figure(figsize=(8, 5))
plt.barh(top10['Family_ID'], top10['Epsilon_Squared'])
plt.xlabel('Kruskal-Wallis Epsilon-Squared')
plt.ylabel('BRCA-PRISM Family')
plt.title('Top PAM50-Associated BRCA-PRISM Families')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'BRCA_PRISM_step9C_top10_effect_sizes.png', dpi=300, bbox_inches='tight')
plt.show()
top_ids = association.sort_values('Epsilon_Squared', ascending=False).head(10)['Family_ID']
heat = association.set_index('Family_ID').loc[top_ids, ['Median_Basal', 'Median_Her2', 'Median_LumA', 'Median_LumB', 'Median_Normal-like']]
heat.columns = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
plt.figure(figsize=(8, 6))
plt.imshow(heat.values, aspect='auto')
plt.xticks(range(len(heat.columns)), heat.columns)
plt.yticks(range(len(heat.index)), heat.index)
plt.colorbar(label='Median standardized family activity')
plt.xlabel('PAM50 Subtype')
plt.ylabel('BRCA-PRISM Family')
plt.title('Subtype-Specific Activity of Top BRCA-PRISM Families')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'BRCA_PRISM_step9C_PAM50_activity_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()
print('STEP 9C COMPLETE')
print('\nTop 10 families:')
display(association[['Family_ID', 'Epsilon_Squared', 'FDR']].sort_values('Epsilon_Squared', ascending=False).head(10).round(5))
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
ASSOC_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
families = pd.read_csv(FAMILY_FILE)
assoc = pd.read_csv(ASSOC_FILE)
top10_ids = assoc.sort_values('Epsilon_Squared', ascending=False).head(10)['Family_ID']
top10_names = families[families['Family_ID'].isin(top10_ids)].groupby('Family_ID').agg(Pathway_Names=('Pathway_Name', lambda x: ' | '.join(sorted(set(x)))), Resources=('Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
top10_names = assoc[['Family_ID', 'Epsilon_Squared', 'FDR']].merge(top10_names, on='Family_ID', how='left').sort_values('Epsilon_Squared', ascending=False).head(10)
display(top10_names)
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
ASSOC_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
assoc = pd.read_csv(ASSOC_FILE)
top10 = assoc.sort_values('Epsilon_Squared', ascending=False).head(10).copy()
print('STEP 9C')
print('\nTop 10 families by effect size:')
display(top10[['Family_ID', 'Epsilon_Squared', 'FDR', 'Median_Basal', 'Median_Her2', 'Median_LumA', 'Median_LumB', 'Median_Normal-like']].round(5))
plot_df = top10.sort_values('Epsilon_Squared')
plt.figure(figsize=(8, 5))
plt.barh(plot_df['Family_ID'], plot_df['Epsilon_Squared'])
plt.xlabel('Epsilon-Squared')
plt.ylabel('BRCA-PRISM Family')
plt.title('Top PAM50-Associated BRCA-PRISM Families')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'BRCA_PRISM_step9C_top10_effect_sizes.png', dpi=300, bbox_inches='tight')
plt.show()
heat = top10[['Median_Basal', 'Median_Her2', 'Median_LumA', 'Median_LumB', 'Median_Normal-like']].copy()
heat.index = top10['Family_ID']
heat.columns = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
plt.figure(figsize=(8, 6))
plt.imshow(heat.values, aspect='auto')
plt.xticks(range(5), heat.columns)
plt.yticks(range(len(heat)), heat.index)
plt.colorbar(label='Median standardized family activity')
plt.xlabel('PAM50 Subtype')
plt.ylabel('BRCA-PRISM Family')
plt.title('Subtype-Specific Activity of Top BRCA-PRISM Families')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'BRCA_PRISM_step9C_PAM50_activity_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()
print('\nSTEP 9C COMPLETE')
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
assoc = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv')
families = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv')
top10_ids = assoc.sort_values('Epsilon_Squared', ascending=False).head(10)['Family_ID']
top10_names = families[families['Family_ID'].isin(top10_ids)].groupby('Family_ID').agg(Resources=('Database', lambda x: '; '.join(sorted(set(x)))), Pathway_Names=('Pathway_Name', lambda x: ' | '.join(sorted(set(x))))).reset_index()
top10_names = assoc[['Family_ID', 'Epsilon_Squared', 'FDR']].merge(top10_names, on='Family_ID').sort_values('Epsilon_Squared', ascending=False).head(10)
display(top10_names)
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
tcga = pd.read_csv(TCGA_FILE)
families = pd.read_csv(FAMILY_FILE)
master = pd.read_csv(MASTER_FILE)
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master['Pathway_ID'] = master['Pathway_ID'].astype(str)
families['Pathway_ID'] = families['Pathway_ID'].astype(str)
family_resource_genes = families[['Family_ID', 'Database', 'Pathway_ID']].merge(master[['Database', 'Pathway_ID', 'Gene_Symbol']], on=['Database', 'Pathway_ID'], how='left').dropna(subset=['Gene_Symbol']).drop_duplicates()
print('Family-resource gene associations:', len(family_resource_genes))
needed_genes = set(family_resource_genes['Gene_Symbol'])
rna_cols = [c for c in tcga.columns if c.lower().endswith('_rna') and c[:-4].upper() in needed_genes]
cnv_cols = [c for c in tcga.columns if c.lower().endswith('_cnv') and c[:-4].upper() in needed_genes]
feature_cols = rna_cols + cnv_cols
print('RNA features used:', len(rna_cols))
print('CNV features used:', len(cnv_cols))
X = tcga[feature_cols].apply(pd.to_numeric, errors='coerce')
means = X.mean(axis=0)
stds = X.std(axis=0, ddof=0).replace(0, np.nan)
Z = (X - means) / stds
scores = tcga[['Cell_line', 'PAM50_Subtype']].copy()
score_info = []
for (family_id, database), g in family_resource_genes.groupby(['Family_ID', 'Database']):
    genes = set(g['Gene_Symbol'])
    fam_rna = [c for c in rna_cols if c[:-4].upper() in genes]
    fam_cnv = [c for c in cnv_cols if c[:-4].upper() in genes]
    rna_score = Z[fam_rna].mean(axis=1) if fam_rna else pd.Series(np.nan, index=Z.index)
    cnv_score = Z[fam_cnv].mean(axis=1) if fam_cnv else pd.Series(np.nan, index=Z.index)
    combined = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1)
    combined = (combined - combined.mean()) / combined.std(ddof=0)
    col = f'{family_id}__{database}'
    scores[col] = combined
    score_info.append({'Family_ID': family_id, 'Database': database, 'Unique_Genes': len(genes), 'RNA_Features': len(fam_rna), 'CNV_Features': len(fam_cnv), 'Score_Column': col})
score_info = pd.DataFrame(score_info)
results = []
for family_id, group in score_info.groupby('Family_ID'):
    rows = group.reset_index(drop=True)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            a = rows.iloc[i]
            b = rows.iloc[j]
            x = scores[a['Score_Column']]
            y = scores[b['Score_Column']]
            valid = x.notna() & y.notna()
            rho, p = spearmanr(x[valid], y[valid])
            mean_abs_diff = np.mean(np.abs(x[valid].values - y[valid].values))
            results.append({'Family_ID': family_id, 'Resource_A': a['Database'], 'Resource_B': b['Database'], 'Genes_A': a['Unique_Genes'], 'Genes_B': b['Unique_Genes'], 'Spearman_Rho': rho, 'P_Value': p, 'Mean_Absolute_Score_Difference': mean_abs_diff, 'N_Patients': int(valid.sum())})
stability = pd.DataFrame(results)
family_stability = stability.groupby('Family_ID').agg(Number_of_Resource_Pairs=('Spearman_Rho', 'count'), Mean_Spearman_Rho=('Spearman_Rho', 'mean'), Median_Spearman_Rho=('Spearman_Rho', 'median'), Min_Spearman_Rho=('Spearman_Rho', 'min'), Max_Spearman_Rho=('Spearman_Rho', 'max'), Mean_Absolute_Score_Difference=('Mean_Absolute_Score_Difference', 'mean')).reset_index().sort_values('Mean_Spearman_Rho', ascending=False)
print('\n' + '=' * 100)
print('STEP 10A CROSS-RESOURCE STABILITY')
print('=' * 100)
print('\nMost stable families:')
display(family_stability.head(10).round(4))
print('\nLeast stable families:')
display(family_stability.tail(10).sort_values('Mean_Spearman_Rho').round(4))
scores.to_csv(RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv', index=False)
score_info.to_csv(RESULTS_DIR / 'BRCA_PRISM_step10A_resource_family_features.csv', index=False)
stability.to_csv(RESULTS_DIR / 'BRCA_PRISM_step10A_pairwise_resource_stability.csv', index=False)
family_stability.to_csv(RESULTS_DIR / 'BRCA_PRISM_step10A_family_stability_summary.csv', index=False)
print('\nSTEP 10A COMPLETE')
print('\nSend me:')
print('1. Most stable families table')
print('2. Least stable families table')
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
STABILITY_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_pairwise_resource_stability.csv'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
stability = pd.read_csv(STABILITY_FILE)
families = pd.read_csv(FAMILY_FILE)
master = pd.read_csv(MASTER_FILE)
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master['Pathway_ID'] = master['Pathway_ID'].astype(str)
families['Pathway_ID'] = families['Pathway_ID'].astype(str)
family_resource = families[['Family_ID', 'Database', 'Pathway_ID']].merge(master[['Database', 'Pathway_ID', 'Gene_Symbol']], on=['Database', 'Pathway_ID'], how='left').dropna(subset=['Gene_Symbol'])
gene_sets = {}
for (family_id, database), g in family_resource.groupby(['Family_ID', 'Database']):
    gene_sets[family_id, database] = set(g['Gene_Symbol'])
jaccards = []
for _, row in stability.iterrows():
    family_id = row['Family_ID']
    a = row['Resource_A']
    b = row['Resource_B']
    genes_a = gene_sets[family_id, a]
    genes_b = gene_sets[family_id, b]
    shared = len(genes_a & genes_b)
    union = len(genes_a | genes_b)
    jaccard = shared / union if union else 0
    jaccards.append({'Family_ID': family_id, 'Resource_A': a, 'Resource_B': b, 'Shared_Genes': shared, 'Gene_Union': union, 'Gene_Jaccard': jaccard})
jaccards = pd.DataFrame(jaccards)
combined = stability.merge(jaccards, on=['Family_ID', 'Resource_A', 'Resource_B'], how='left')
rho, p = spearmanr(combined['Gene_Jaccard'], combined['Spearman_Rho'])
print('=' * 90)
print('STEP 10B')
print('=' * 90)
print('\nNumber of resource pairs:', len(combined))
print('Spearman correlation between Gene Jaccard and', 'Patient-score stability:')
print('rho =', round(rho, 4))
print('p   =', p)
print('\nLowest-stability resource pairs:')
display(combined[['Family_ID', 'Resource_A', 'Resource_B', 'Gene_Jaccard', 'Spearman_Rho', 'Mean_Absolute_Score_Difference']].sort_values('Spearman_Rho').head(15).round(4))
print('\nHighest-stability resource pairs:')
display(combined[['Family_ID', 'Resource_A', 'Resource_B', 'Gene_Jaccard', 'Spearman_Rho', 'Mean_Absolute_Score_Difference']].sort_values('Spearman_Rho', ascending=False).head(15).round(4))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step10B_gene_overlap_vs_patient_stability.csv'
combined.to_csv(OUT_FILE, index=False)
print('\nSTEP 10B COMPLETE')
print('Saved:', OUT_FILE)
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
FILE = RESULTS_DIR / 'BRCA_PRISM_step10B_gene_overlap_vs_patient_stability.csv'
df = pd.read_csv(FILE)
plt.figure(figsize=(7, 6))
plt.scatter(df['Gene_Jaccard'], df['Spearman_Rho'], s=55, alpha=0.8)
plt.xlabel('Gene-set Jaccard similarity')
plt.ylabel('Patient-level Spearman correlation')
plt.title('Gene-Set Similarity vs Patient-Level Stability')
plt.text(0.05, 0.95, '$\\rho = 0.851$\n$p = 5.27 \\times 10^{-17}$', transform=plt.gca().transAxes, verticalalignment='top')
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'BRCA_PRISM_step10C_gene_overlap_vs_patient_stability.png', dpi=300, bbox_inches='tight')
plt.show()
print('STEP 10C COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kruskal
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
SCORES_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
INFO_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_family_features.csv'
scores = pd.read_csv(SCORES_FILE)
info = pd.read_csv(INFO_FILE)
subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
results = []
for _, row in info.iterrows():
    family_id = row['Family_ID']
    database = row['Database']
    col = row['Score_Column']
    groups = [scores.loc[scores['PAM50_Subtype'] == subtype, col].dropna().values for subtype in subtypes]
    H, p = kruskal(*groups)
    n = sum((len(g) for g in groups))
    k = len(groups)
    epsilon2 = max(0, (H - k + 1) / (n - k))
    results.append({'Family_ID': family_id, 'Database': database, 'Kruskal_H': H, 'P_Value': p, 'Epsilon_Squared': epsilon2})
resource_assoc = pd.DataFrame(results)
resource_assoc = resource_assoc.sort_values('P_Value').reset_index(drop=True)
m = len(resource_assoc)
resource_assoc['FDR'] = resource_assoc['P_Value'] * m / np.arange(1, m + 1)
resource_assoc['FDR'] = resource_assoc['FDR'][::-1].cummin()[::-1].clip(upper=1)
effect_stability = resource_assoc.groupby('Family_ID').agg(Number_of_Resources=('Database', 'nunique'), Mean_Epsilon_Squared=('Epsilon_Squared', 'mean'), Min_Epsilon_Squared=('Epsilon_Squared', 'min'), Max_Epsilon_Squared=('Epsilon_Squared', 'max'), SD_Epsilon_Squared=('Epsilon_Squared', 'std')).reset_index()
effect_stability['Effect_Size_Range'] = effect_stability['Max_Epsilon_Squared'] - effect_stability['Min_Epsilon_Squared']
effect_stability = effect_stability.sort_values('Effect_Size_Range', ascending=False)
print('=' * 100)
print('STEP 11 RESOURCE-SPECIFIC PAM50 DISCRIMINATION')
print('=' * 100)
print('\nFamilies with largest resource-dependent effect-size variation:')
display(effect_stability.head(15).round(4))
print('\nTop resource-specific PAM50 associations:')
display(resource_assoc[['Family_ID', 'Database', 'Epsilon_Squared', 'FDR']].sort_values('Epsilon_Squared', ascending=False).head(20).round(5))
resource_assoc.to_csv(RESULTS_DIR / 'BRCA_PRISM_step11_resource_specific_PAM50_association.csv', index=False)
effect_stability.to_csv(RESULTS_DIR / 'BRCA_PRISM_step11_family_effect_size_stability.csv', index=False)
print('\nSTEP 11 COMPLETE')
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step11_resource_specific_PAM50_association.csv'
consensus = pd.read_csv(CONSENSUS_FILE)
resource = pd.read_csv(RESOURCE_FILE)
resource_summary = resource.groupby('Family_ID').agg(Number_of_Resources=('Database', 'nunique'), Mean_Resource_Effect=('Epsilon_Squared', 'mean'), Min_Resource_Effect=('Epsilon_Squared', 'min'), Max_Resource_Effect=('Epsilon_Squared', 'max'), SD_Resource_Effect=('Epsilon_Squared', 'std')).reset_index()
comparison = consensus[['Family_ID', 'Epsilon_Squared']].rename(columns={'Epsilon_Squared': 'Consensus_Effect'}).merge(resource_summary, on='Family_ID', how='left')
comparison['Consensus_minus_Mean'] = comparison['Consensus_Effect'] - comparison['Mean_Resource_Effect']
comparison['Consensus_minus_Max'] = comparison['Consensus_Effect'] - comparison['Max_Resource_Effect']
comparison['Consensus_within_resource_range'] = (comparison['Consensus_Effect'] >= comparison['Min_Resource_Effect']) & (comparison['Consensus_Effect'] <= comparison['Max_Resource_Effect'])
print('=' * 100)
print('STEP 12 CONSENSUS VS RESOURCE-SPECIFIC PAM50 EFFECT')
print('=' * 100)
print('\nConsensus effect within individual-resource range:', comparison['Consensus_within_resource_range'].sum(), '/', len(comparison))
print('Consensus effect >= mean resource effect:', (comparison['Consensus_Effect'] >= comparison['Mean_Resource_Effect']).sum(), '/', len(comparison))
print('\nLargest consensus advantages over mean resource effect:')
display(comparison[['Family_ID', 'Consensus_Effect', 'Mean_Resource_Effect', 'Min_Resource_Effect', 'Max_Resource_Effect', 'Consensus_minus_Mean']].sort_values('Consensus_minus_Mean', ascending=False).head(15).round(4))
print('\nLargest consensus disadvantages:')
display(comparison[['Family_ID', 'Consensus_Effect', 'Mean_Resource_Effect', 'Min_Resource_Effect', 'Max_Resource_Effect', 'Consensus_minus_Mean']].sort_values('Consensus_minus_Mean').head(15).round(4))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step12_consensus_vs_resource_effects.csv'
comparison.to_csv(OUT_FILE, index=False)
print('\nSTEP 12 COMPLETE')
print('Saved:', OUT_FILE)
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_SCORE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_SCORE_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_SCORE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_SCORE_FILE)
print('Resource score table:', resource_scores.shape)
print('Consensus score table:', consensus_scores.shape)
consensus_cols = [c for c in consensus_scores.columns if c.endswith('_Combined')]
print('Consensus family features:', len(consensus_cols))
resources = ['KEGG', 'Reactome', 'WikiPathways', 'LCPathways']
datasets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    if len(cols) > 0:
        datasets[resource] = resource_scores[cols].copy()
        print(f'{resource}:', len(cols), 'family features')
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy()
y = resource_scores['PAM50_Subtype'].copy()
print('\nSubtype counts:')
print(y.value_counts())
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
all_results = []
for representation, X in datasets.items():
    X = X.dropna(axis=1, how='all')
    X = X.fillna(X.mean())
    print(f'\nRunning {representation}:', X.shape[1], 'features')
    fold_number = 1
    for train_idx, test_idx in cv.split(X, y):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        all_results.append({'Representation': representation, 'Fold': fold_number, 'Number_of_Features': X.shape[1], 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
        fold_number += 1
cv_results = pd.DataFrame(all_results)
summary = cv_results.groupby('Representation').agg(Number_of_Features=('Number_of_Features', 'first'), Mean_Accuracy=('Accuracy', 'mean'), SD_Accuracy=('Accuracy', 'std'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), SD_Balanced_Accuracy=('Balanced_Accuracy', 'std'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index()
summary = summary.sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13 PREDICTION STABILITY')
print('=' * 100)
display(summary.round(4))
print('\nFold-level results:')
display(cv_results.round(4))
cv_results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13_prediction_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13_prediction_summary.csv', index=False)
print('\nSTEP 13 COMPLETE')
print('\nSend me:')
print('1. STEP 13 PREDICTION STABILITY summary table')
from pathlib import Path
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways', 'LCPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
    print(resource, 'families:', len(family_sets[resource]))
common_families = set.intersection(*family_sets.values())
common_families = sorted(common_families)
print('\nFamilies common to all 4 resources:')
print(common_families)
print('\nNumber of matched families:', len(common_families))
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in common_families]
    datasets[resource] = resource_scores[cols].copy()
consensus_cols = [f'{family}_Combined' for family in common_families]
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy()
y = resource_scores['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    X = X.fillna(X.mean())
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Representation': representation, 'Fold': fold, 'Matched_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Matched_Families=('Matched_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13B MATCHED-FAMILY COMPARISON')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13B_matched_family_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13B_matched_family_summary.csv', index=False)
print('\nSTEP 13B COMPLETE')
from pathlib import Path
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
    print(resource, 'families:', len(family_sets[resource]))
common_families = sorted(set.intersection(*family_sets.values()))
print('\nCommon families:')
print(common_families)
print('\nNumber of matched families:', len(common_families))
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in common_families]
    datasets[resource] = resource_scores[cols].copy()
consensus_cols = [f'{family}_Combined' for family in common_families]
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy()
y = resource_scores['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    X = X.fillna(X.mean())
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Representation': representation, 'Fold': fold, 'Matched_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Matched_Families=('Matched_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13C MATCHED 3-RESOURCE COMPARISON')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13C_matched_3resource_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13C_matched_3resource_summary.csv', index=False)
print('\nSTEP 13C COMPLETE')
from pathlib import Path
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
common_families = sorted(set.intersection(*family_sets.values()))
print('Matched families:', len(common_families))
print(common_families)
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in common_families]
    datasets[resource] = resource_scores[cols].copy().fillna(0)
consensus_cols = [f'{family}_Combined' for family in common_families]
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy().fillna(0)
y = resource_scores['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    print('\nRunning:', representation)
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = SVC(kernel='linear', class_weight='balanced', C=1.0)
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Model': 'Linear_SVM', 'Representation': representation, 'Fold': fold, 'Matched_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Matched_Families=('Matched_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13D LINEAR SVM')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13D_linearSVM_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13D_linearSVM_summary.csv', index=False)
print('\nSTEP 13D COMPLETE')
from pathlib import Path
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
common_families = sorted(set.intersection(*family_sets.values()))
print('Matched families:', len(common_families))
print(common_families)
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in common_families]
    datasets[resource] = resource_scores[cols].copy().fillna(0)
consensus_cols = [f'{family}_Combined' for family in common_families]
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy().fillna(0)
y = resource_scores['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    print('\nRunning:', representation)
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = RandomForestClassifier(n_estimators=500, class_weight='balanced', random_state=42, n_jobs=-1)
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Model': 'Random_Forest', 'Representation': representation, 'Fold': fold, 'Matched_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Matched_Families=('Matched_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13E RANDOM FOREST')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13E_randomForest_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13E_randomForest_summary.csv', index=False)
print('\nSTEP 13E COMPLETE')
from pathlib import Path
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
common_families = sorted(set.intersection(*family_sets.values()))
print('Matched families:', len(common_families))
print(common_families)
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in common_families]
    datasets[resource] = resource_scores[cols].copy().fillna(0)
consensus_cols = [f'{family}_Combined' for family in common_families]
datasets['BRCA_PRISM_Consensus'] = consensus_scores[consensus_cols].copy().fillna(0)
y = resource_scores['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    print('\nRunning:', representation)
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_depth=3, random_state=42)
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Model': 'Gradient_Boosting', 'Representation': representation, 'Fold': fold, 'Matched_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Matched_Families=('Matched_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 13F GRADIENT BOOSTING')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13F_gradientBoosting_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step13F_gradientBoosting_summary.csv', index=False)
print('\nSTEP 13F COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
matched_families = sorted(set.intersection(*family_sets.values()))
print('Matched families:', len(matched_families))
print(matched_families)
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in matched_families]
    X = resource_scores[cols].copy()
    X.columns = matched_families
    datasets[resource] = X
consensus_cols = [f'{family}_Combined' for family in matched_families]
X_consensus = consensus_scores[consensus_cols].copy()
X_consensus.columns = matched_families
datasets['BRCA_PRISM_Consensus'] = X_consensus
y = resource_scores['PAM50_Subtype'].copy()
importance_rows = []
for representation, X in datasets.items():
    X = X.fillna(X.mean())
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
    model.fit(X_scaled, y)
    importance = np.sqrt(np.sum(model.coef_ ** 2, axis=0))
    temp = pd.DataFrame({'Representation': representation, 'Family_ID': matched_families, 'Importance': importance})
    temp['Rank'] = temp['Importance'].rank(ascending=False, method='min').astype(int)
    importance_rows.append(temp)
importance_table = pd.concat(importance_rows, ignore_index=True)
print('\n' + '=' * 100)
print('STEP 14A PATHWAY IMPORTANCE RANKINGS')
print('=' * 100)
for representation in datasets.keys():
    print('\n', representation)
    display(importance_table[importance_table['Representation'] == representation].sort_values('Rank').round(4))
rank_matrix = importance_table.pivot(index='Family_ID', columns='Representation', values='Rank')
print('\nRank matrix:')
display(rank_matrix)
representations = list(rank_matrix.columns)
agreement_rows = []
for i in range(len(representations)):
    for j in range(i + 1, len(representations)):
        a = representations[i]
        b = representations[j]
        rho, p = spearmanr(rank_matrix[a], rank_matrix[b])
        agreement_rows.append({'Representation_A': a, 'Representation_B': b, 'Spearman_Rank_Rho': rho, 'P_Value': p})
rank_agreement = pd.DataFrame(agreement_rows)
print('\nPairwise pathway-ranking agreement:')
display(rank_agreement.sort_values('Spearman_Rank_Rho', ascending=False).round(4))
top3 = {}
for representation in datasets.keys():
    top3[representation] = set(importance_table[importance_table['Representation'] == representation].sort_values('Rank').head(3)['Family_ID'])
top3_rows = []
for i in range(len(representations)):
    for j in range(i + 1, len(representations)):
        a = representations[i]
        b = representations[j]
        shared = len(top3[a] & top3[b])
        union = len(top3[a] | top3[b])
        top3_rows.append({'Representation_A': a, 'Representation_B': b, 'Shared_Top3_Families': shared, 'Top3_Jaccard': shared / union if union > 0 else 0})
top3_agreement = pd.DataFrame(top3_rows)
print('\nTop-3 family overlap:')
display(top3_agreement.sort_values('Top3_Jaccard', ascending=False).round(4))
importance_table.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14A_logistic_family_importance.csv', index=False)
rank_agreement.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14A_logistic_rank_agreement.csv', index=False)
top3_agreement.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14A_logistic_top3_overlap.csv', index=False)
print('\nSTEP 14A COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
RESOURCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step10A_resource_specific_family_scores.csv'
CONSENSUS_FILE = RESULTS_DIR / 'BRCA_PRISM_step9B_TCGA_family_activity_scores.csv'
resource_scores = pd.read_csv(RESOURCE_FILE)
consensus_scores = pd.read_csv(CONSENSUS_FILE)
resources = ['KEGG', 'Reactome', 'WikiPathways']
family_sets = {}
for resource in resources:
    cols = [c for c in resource_scores.columns if c.endswith(f'__{resource}')]
    family_sets[resource] = {c.split('__')[0] for c in cols}
matched_families = sorted(set.intersection(*family_sets.values()))
print('Matched families:', len(matched_families))
print(matched_families)
datasets = {}
for resource in resources:
    cols = [f'{family}__{resource}' for family in matched_families]
    X = resource_scores[cols].copy()
    X.columns = matched_families
    datasets[resource] = X
consensus_cols = [f'{family}_Combined' for family in matched_families]
X_consensus = consensus_scores[consensus_cols].copy()
X_consensus.columns = matched_families
datasets['BRCA_PRISM_Consensus'] = X_consensus
y = resource_scores['PAM50_Subtype'].copy()
sample_ids = resource_scores['Cell_line'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
splits = list(cv.split(datasets['KEGG'], y))
explanation_rows = []
for representation, X in datasets.items():
    X = X.copy()
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx].copy()
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        train_means = X_train.mean()
        X_train = X_train.fillna(train_means)
        X_test = X_test.fillna(train_means)
        scaler = StandardScaler()
        X_train_z = scaler.fit_transform(X_train)
        X_test_z = scaler.transform(X_test)
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(X_train_z, y_train)
        predictions = model.predict(X_test_z)
        class_to_index = {cls: i for i, cls in enumerate(model.classes_)}
        for local_i, sample_index in enumerate(test_idx):
            predicted_class = predictions[local_i]
            coef_row = model.coef_[class_to_index[predicted_class]]
            contributions = np.abs(X_test_z[local_i] * coef_row)
            order = np.argsort(-contributions)
            ranked_families = [matched_families[i] for i in order]
            row = {'Cell_line': sample_ids.iloc[sample_index], 'True_Subtype': y_test.loc[y_test.index[local_i]], 'Representation': representation, 'Predicted_Subtype': predicted_class, 'Fold': fold, 'Top1_Family': ranked_families[0], 'Top2_Family': ranked_families[1], 'Top3_Family': ranked_families[2]}
            for family, value in zip(matched_families, contributions):
                row[f'Contribution__{family}'] = value
            explanation_rows.append(row)
explanations = pd.DataFrame(explanation_rows)
print('\nPatient-explanation rows:', len(explanations))
representations = list(datasets.keys())
pair_rows = []
for i in range(len(representations)):
    for j in range(i + 1, len(representations)):
        rep_a = representations[i]
        rep_b = representations[j]
        A = explanations[explanations['Representation'] == rep_a].set_index('Cell_line')
        B = explanations[explanations['Representation'] == rep_b].set_index('Cell_line')
        shared_ids = A.index.intersection(B.index)
        for sid in shared_ids:
            top_a = {A.loc[sid, 'Top1_Family'], A.loc[sid, 'Top2_Family'], A.loc[sid, 'Top3_Family']}
            top_b = {B.loc[sid, 'Top1_Family'], B.loc[sid, 'Top2_Family'], B.loc[sid, 'Top3_Family']}
            shared_top3 = len(top_a & top_b)
            top3_jaccard = shared_top3 / len(top_a | top_b)
            contrib_a = np.array([A.loc[sid, f'Contribution__{family}'] for family in matched_families])
            contrib_b = np.array([B.loc[sid, f'Contribution__{family}'] for family in matched_families])
            rho, _ = spearmanr(contrib_a, contrib_b)
            pair_rows.append({'Cell_line': sid, 'Representation_A': rep_a, 'Representation_B': rep_b, 'Prediction_A': A.loc[sid, 'Predicted_Subtype'], 'Prediction_B': B.loc[sid, 'Predicted_Subtype'], 'Same_Predicted_Subtype': A.loc[sid, 'Predicted_Subtype'] == B.loc[sid, 'Predicted_Subtype'], 'Shared_Top3_Families': shared_top3, 'Top3_Jaccard': top3_jaccard, 'Contribution_Rank_Rho': rho})
patient_agreement = pd.DataFrame(pair_rows)
summary = patient_agreement.groupby(['Representation_A', 'Representation_B']).agg(N_Patients=('Cell_line', 'count'), Prediction_Agreement=('Same_Predicted_Subtype', 'mean'), Mean_Shared_Top3=('Shared_Top3_Families', 'mean'), Mean_Top3_Jaccard=('Top3_Jaccard', 'mean'), Mean_Contribution_Rank_Rho=('Contribution_Rank_Rho', 'mean')).reset_index()
same_prediction_summary = patient_agreement[patient_agreement['Same_Predicted_Subtype']].groupby(['Representation_A', 'Representation_B']).agg(N_Same_Prediction=('Cell_line', 'count'), Mean_Shared_Top3=('Shared_Top3_Families', 'mean'), Mean_Top3_Jaccard=('Top3_Jaccard', 'mean'), Mean_Contribution_Rank_Rho=('Contribution_Rank_Rho', 'mean')).reset_index()
print('\n' + '=' * 100)
print('STEP 14B PATIENT-LEVEL EXPLANATION STABILITY')
print('=' * 100)
print('\nAll patients:')
display(summary.round(4))
print('\nAmong patients where BOTH representations predicted the SAME PAM50 subtype:')
display(same_prediction_summary.round(4))
identical_top3 = patient_agreement['Shared_Top3_Families'] == 3
print('\nOverall fraction with identical Top-3 families:', round(identical_top3.mean(), 4))
explanations.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14B_patient_explanations.csv', index=False)
patient_agreement.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14B_patient_explanation_pairwise_agreement.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14B_patient_explanation_summary.csv', index=False)
same_prediction_summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14B_same_prediction_explanation_summary.csv', index=False)
print('\nSTEP 14B COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
IMPORTANCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step14A_logistic_family_importance.csv'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
CONSENSUS_GENE_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
importance = pd.read_csv(IMPORTANCE_FILE)
families = pd.read_csv(FAMILY_FILE)
master = pd.read_csv(MASTER_FILE)
consensus_genes = pd.read_csv(CONSENSUS_GENE_FILE)
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master['Pathway_ID'] = master['Pathway_ID'].astype(str)
families['Pathway_ID'] = families['Pathway_ID'].astype(str)
consensus_genes['Gene_Symbol'] = consensus_genes['Gene_Symbol'].astype(str).str.strip().str.upper()
resource_family_genes = families[['Family_ID', 'Database', 'Pathway_ID']].merge(master[['Database', 'Pathway_ID', 'Gene_Symbol']], on=['Database', 'Pathway_ID'], how='left').dropna(subset=['Gene_Symbol']).drop_duplicates(subset=['Family_ID', 'Database', 'Gene_Symbol'])
gene_importance_rows = []
representations = ['KEGG', 'Reactome', 'WikiPathways', 'BRCA_PRISM_Consensus']
for representation in representations:
    imp = importance[importance['Representation'] == representation].copy()
    for _, row in imp.iterrows():
        family_id = row['Family_ID']
        family_importance = row['Importance']
        if representation == 'BRCA_PRISM_Consensus':
            genes = consensus_genes[(consensus_genes['Family_ID'] == family_id) & (consensus_genes['Resource_Support'] >= 2)]['Gene_Symbol'].dropna().unique().tolist()
        else:
            genes = resource_family_genes[(resource_family_genes['Family_ID'] == family_id) & (resource_family_genes['Database'] == representation)]['Gene_Symbol'].dropna().unique().tolist()
        if len(genes) == 0:
            continue
        per_gene_weight = family_importance / len(genes)
        for gene in genes:
            gene_importance_rows.append({'Representation': representation, 'Family_ID': family_id, 'Gene_Symbol': gene, 'Family_Importance': family_importance, 'Family_Gene_Count': len(genes), 'Gene_Contribution': per_gene_weight})
gene_long = pd.DataFrame(gene_importance_rows)
gene_rank = gene_long.groupby(['Representation', 'Gene_Symbol']).agg(Gene_Importance=('Gene_Contribution', 'sum'), Number_of_Important_Families=('Family_ID', 'nunique')).reset_index()
gene_rank['Gene_Rank'] = gene_rank.groupby('Representation')['Gene_Importance'].rank(ascending=False, method='min')
print('\n' + '=' * 100)
print('STEP 14C GENE-LEVEL INTERPRETATION')
print('=' * 100)
for representation in representations:
    print('\n', representation)
    display(gene_rank[gene_rank['Representation'] == representation].sort_values('Gene_Rank').head(20).round(6))
top20 = {}
for representation in representations:
    top20[representation] = set(gene_rank[gene_rank['Representation'] == representation].sort_values('Gene_Rank').head(20)['Gene_Symbol'])
overlap_rows = []
for i in range(len(representations)):
    for j in range(i + 1, len(representations)):
        a = representations[i]
        b = representations[j]
        shared = len(top20[a] & top20[b])
        union = len(top20[a] | top20[b])
        overlap_rows.append({'Representation_A': a, 'Representation_B': b, 'Shared_Top20_Genes': shared, 'Top20_Jaccard': shared / union if union else 0})
top20_overlap = pd.DataFrame(overlap_rows)
print('\nTop-20 gene overlap:')
display(top20_overlap.sort_values('Top20_Jaccard', ascending=False).round(4))
rank_rows = []
for i in range(len(representations)):
    for j in range(i + 1, len(representations)):
        a = representations[i]
        b = representations[j]
        A = gene_rank[gene_rank['Representation'] == a][['Gene_Symbol', 'Gene_Importance']].rename(columns={'Gene_Importance': 'Importance_A'})
        B = gene_rank[gene_rank['Representation'] == b][['Gene_Symbol', 'Gene_Importance']].rename(columns={'Gene_Importance': 'Importance_B'})
        shared = A.merge(B, on='Gene_Symbol', how='inner')
        if len(shared) >= 3:
            rho, p = spearmanr(shared['Importance_A'], shared['Importance_B'])
        else:
            rho = np.nan
            p = np.nan
        rank_rows.append({'Representation_A': a, 'Representation_B': b, 'Shared_Genes': len(shared), 'Gene_Importance_Rho': rho, 'P_Value': p})
gene_rank_agreement = pd.DataFrame(rank_rows)
print('\nGene-importance ranking agreement:')
display(gene_rank_agreement.sort_values('Gene_Importance_Rho', ascending=False).round(4))
top20_membership = []
all_top_genes = sorted(set.union(*top20.values()))
for gene in all_top_genes:
    reps = [rep for rep in representations if gene in top20[rep]]
    top20_membership.append({'Gene_Symbol': gene, 'Number_of_Representations': len(reps), 'Representations': '; '.join(reps)})
consensus_top_genes = pd.DataFrame(top20_membership).sort_values(['Number_of_Representations', 'Gene_Symbol'], ascending=[False, True])
print('\nCross-resource consensus Top-20 genes:')
display(consensus_top_genes.head(30))
gene_rank.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14C_gene_importance_rankings.csv', index=False)
top20_overlap.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14C_top20_gene_overlap.csv', index=False)
gene_rank_agreement.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14C_gene_rank_agreement.csv', index=False)
consensus_top_genes.to_csv(RESULTS_DIR / 'BRCA_PRISM_step14C_consensus_top_genes.csv', index=False)
print('\nSTEP 14C COMPLETE')
from pathlib import Path
import pandas as pd
import numpy as np
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
IMPORTANCE_FILE = RESULTS_DIR / 'BRCA_PRISM_step14A_logistic_family_importance.csv'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
importance = pd.read_csv(IMPORTANCE_FILE)
families = pd.read_csv(FAMILY_FILE)
family_names = families.groupby('Family_ID').agg(Pathway_Names=('Pathway_Name', lambda x: ' | '.join(sorted(set(x)))), Resources_in_Family=('Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
summary = importance.groupby('Family_ID').agg(Mean_Rank=('Rank', 'mean'), Median_Rank=('Rank', 'median'), Best_Rank=('Rank', 'min'), Worst_Rank=('Rank', 'max'), Rank_SD=('Rank', 'std'), Mean_Importance=('Importance', 'mean'), Max_Importance=('Importance', 'max'), Min_Importance=('Importance', 'min')).reset_index()
summary['Rank_Range'] = summary['Worst_Rank'] - summary['Best_Rank']
top3_count = importance[importance['Rank'] <= 3].groupby('Family_ID')['Representation'].nunique().reset_index(name='Top3_Representation_Count')
summary = summary.merge(top3_count, on='Family_ID', how='left')
summary['Top3_Representation_Count'] = summary['Top3_Representation_Count'].fillna(0).astype(int)

def classify(row):
    if row['Top3_Representation_Count'] >= 3:
        return 'Consensus-high importance'
    elif row['Top3_Representation_Count'] == 1:
        return 'Resource-specific high importance'
    else:
        return 'Intermediate / mixed importance'
summary['Importance_Category'] = summary.apply(classify, axis=1)
top3_rows = importance[importance['Rank'] <= 3][['Family_ID', 'Representation', 'Rank']].copy()
top3_resources = top3_rows.groupby('Family_ID').agg(Top3_Representations=('Representation', lambda x: '; '.join(sorted(set(x))))).reset_index()
summary = summary.merge(top3_resources, on='Family_ID', how='left')
summary['Top3_Representations'] = summary['Top3_Representations'].fillna('None')
summary = summary.merge(family_names, on='Family_ID', how='left')
summary = summary.sort_values(['Top3_Representation_Count', 'Mean_Rank'], ascending=[False, True])
print('=' * 100)
print('STEP 14D CONSENSUS VS RESOURCE-SPECIFIC PATHWAY IMPORTANCE')
print('=' * 100)
print('\nConsensus-high importance pathways:')
display(summary[summary['Importance_Category'] == 'Consensus-high importance'][['Family_ID', 'Mean_Rank', 'Rank_Range', 'Top3_Representation_Count', 'Top3_Representations', 'Pathway_Names']].round(3))
print('\nResource-specific high-importance pathways:')
display(summary[summary['Importance_Category'] == 'Resource-specific high importance'][['Family_ID', 'Mean_Rank', 'Rank_Range', 'Top3_Representation_Count', 'Top3_Representations', 'Pathway_Names']].round(3))
print('\nFull pathway importance stability table:')
display(summary[['Family_ID', 'Mean_Rank', 'Best_Rank', 'Worst_Rank', 'Rank_Range', 'Top3_Representation_Count', 'Importance_Category']].round(3))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step14D_consensus_resource_specific_pathway_importance.csv'
summary.to_csv(OUT_FILE, index=False)
print('\nSTEP 14D COMPLETE')
print('Saved:', OUT_FILE)
from pathlib import Path
import pandas as pd
import numpy as np
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
gene_patient = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step10B_gene_overlap_vs_patient_stability.csv')
effect = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step11_family_effect_size_stability.csv')
rank = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step14D_consensus_resource_specific_pathway_importance.csv')
gp_summary = gene_patient.groupby('Family_ID').agg(Mean_Gene_Jaccard=('Gene_Jaccard', 'mean'), Mean_Patient_Spearman=('Spearman_Rho', 'mean')).reset_index()
matched_families = sorted(rank['Family_ID'].unique())
print('Matched families:', len(matched_families))
print(matched_families)
score = rank[['Family_ID', 'Mean_Rank', 'Rank_Range', 'Importance_Category', 'Pathway_Names']].merge(gp_summary, on='Family_ID', how='left').merge(effect[['Family_ID', 'Effect_Size_Range']], on='Family_ID', how='left')
score['Gene_Stability'] = score['Mean_Gene_Jaccard'].clip(0, 1)
score['Patient_Stability'] = score['Mean_Patient_Spearman'].clip(0, 1)
score['Subtype_Effect_Stability'] = (1 - score['Effect_Size_Range']).clip(0, 1)
MAX_RANK_RANGE = len(matched_families) - 1
score['Pathway_Rank_Stability'] = (1 - score['Rank_Range'] / MAX_RANK_RANGE).clip(0, 1)
components = ['Gene_Stability', 'Patient_Stability', 'Subtype_Effect_Stability', 'Pathway_Rank_Stability']
score['BRCA_PRISM_Stability_Score'] = score[components].mean(axis=1)

def stability_label(x):
    if x >= 0.8:
        return 'High'
    elif x >= 0.6:
        return 'Moderate'
    else:
        return 'Low'
score['Descriptive_Stability'] = score['BRCA_PRISM_Stability_Score'].apply(stability_label)
score = score.sort_values('BRCA_PRISM_Stability_Score', ascending=False).reset_index(drop=True)
print('\n' + '=' * 100)
print('STEP 15 BRCA-PRISM STABILITY SCORE')
print('=' * 100)
display(score[['Family_ID', 'Gene_Stability', 'Patient_Stability', 'Subtype_Effect_Stability', 'Pathway_Rank_Stability', 'BRCA_PRISM_Stability_Score', 'Descriptive_Stability']].round(4))
print('\nMost stable families:')
display(score[['Family_ID', 'BRCA_PRISM_Stability_Score', 'Pathway_Names']].head(5).round(4))
print('\nLeast stable families:')
display(score[['Family_ID', 'BRCA_PRISM_Stability_Score', 'Pathway_Names']].tail(5).sort_values('BRCA_PRISM_Stability_Score').round(4))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step15_exploratory_stability_score.csv'
score.to_csv(OUT_FILE, index=False)
print('\nSTEP 15 COMPLETE')
print('Saved:', OUT_FILE)
