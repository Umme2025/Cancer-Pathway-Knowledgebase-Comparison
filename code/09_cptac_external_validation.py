from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
CPTAC_DIR = PROJECT_ROOT / 'data_raw' / 'CPTAC'
FEATURE_FILE = CPTAC_DIR / 'cptac_as_validation.csv'
LABEL_FILE = CPTAC_DIR / 'cptac_brca_cnv_rna_subtypes_independent.csv'
print('Feature file exists:', FEATURE_FILE.exists())
print('Label file exists:', LABEL_FILE.exists())
cptac_x = pd.read_csv(FEATURE_FILE)
cptac_y = pd.read_csv(LABEL_FILE)
print('\nCPTAC feature shape:')
print(cptac_x.shape)
print('\nCPTAC label shape:')
print(cptac_y.shape)
print('\nFeature file first 10 columns:')
print(cptac_x.columns[:10].tolist())
print('\nFeature file last 10 columns:')
print(cptac_x.columns[-10:].tolist())
print('\nLabel file columns:')
print(cptac_y.columns.tolist())
print('\nFirst 5 label rows:')
display(cptac_y.head())
rna_cols = [c for c in cptac_x.columns if c.lower().endswith('_rna')]
cnv_cols = [c for c in cptac_x.columns if c.lower().endswith('_cnv')]
print('\nRNA features:', len(rna_cols))
print('CNV features:', len(cnv_cols))
print('\nExample RNA columns:')
print(rna_cols[:10])
print('\nExample CNV columns:')
print(cnv_cols[:10])
print('\nMissingness:')
print(cptac_x[rna_cols + cnv_cols].isna().mean().mean())
import pandas as pd
print('Feature samples:', cptac_x['Cell_line'].nunique())
print('Label samples:', cptac_y['Cell_line'].nunique())
feature_ids = set(cptac_x['Cell_line'])
label_ids = set(cptac_y['Cell_line'])
print('\nSamples present in both files:', len(feature_ids & label_ids))
print('Feature samples without label:', len(feature_ids - label_ids))
print('Label samples without features:', len(label_ids - feature_ids))
cptac_full = cptac_x.merge(cptac_y, on='Cell_line', how='inner')
print('\nMerged CPTAC shape:')
print(cptac_full.shape)
print('\nPAM50 subtype distribution:')
print(cptac_full['Cancer_type'].value_counts())
print('\nPAM50 subtype percentages:')
print(cptac_full['Cancer_type'].value_counts(normalize=True).mul(100).round(2))
print('\nMissing PAM50 labels:', cptac_full['Cancer_type'].isna().sum())
display(cptac_full[['Cell_line', 'Cancer_type']].head(10))
print(family_gene.columns.tolist())
display(family_gene.head())
import pandas as pd
import numpy as np
cptac_rna_genes = {c[:-4] for c in cptac_x.columns if c.endswith('_rna')}
cptac_cnv_genes = {c[:-4] for c in cptac_x.columns if c.endswith('_cnv')}
cptac_any_genes = cptac_rna_genes | cptac_cnv_genes
print('CPTAC RNA genes:', len(cptac_rna_genes))
print('CPTAC CNV genes:', len(cptac_cnv_genes))
print('CPTAC genes in RNA or CNV:', len(cptac_any_genes))
print('\nBRCA-PRISM families:', family_gene['Family_ID'].nunique())
family_gene['CPTAC_RNA'] = family_gene['Gene_Symbol'].isin(cptac_rna_genes)
family_gene['CPTAC_CNV'] = family_gene['Gene_Symbol'].isin(cptac_cnv_genes)
family_gene['CPTAC_Any'] = family_gene['Gene_Symbol'].isin(cptac_any_genes)
rows = []
for family_id, df in family_gene.groupby('Family_ID'):
    union_genes = set(df['Gene_Symbol'])
    consensus_genes = set(df.loc[df['Resource_Support'] >= 2, 'Gene_Symbol'])
    core_genes = set(df.loc[df['Gene_Category'] == 'Core', 'Gene_Symbol'])
    specific_genes = set(df.loc[df['Gene_Category'] == 'Resource-specific', 'Gene_Symbol'])

    def calc_coverage(genes):
        if len(genes) == 0:
            return (np.nan, 0)
        measured = len(genes.intersection(cptac_any_genes))
        return (measured / len(genes), measured)
    union_cov, union_measured = calc_coverage(union_genes)
    consensus_cov, consensus_measured = calc_coverage(consensus_genes)
    core_cov, core_measured = calc_coverage(core_genes)
    specific_cov, specific_measured = calc_coverage(specific_genes)
    rows.append({'Family_ID': family_id, 'Union_Genes': len(union_genes), 'Union_Measured': union_measured, 'Union_Coverage': union_cov, 'Consensus_Genes': len(consensus_genes), 'Consensus_Measured': consensus_measured, 'Consensus_Coverage': consensus_cov, 'Core_Genes': len(core_genes), 'Core_Measured': core_measured, 'Core_Coverage': core_cov, 'ResourceSpecific_Genes': len(specific_genes), 'ResourceSpecific_Measured': specific_measured, 'ResourceSpecific_Coverage': specific_cov})
coverage_df = pd.DataFrame(rows)
print('\n======================================')
print('CPTAC FAMILY COVERAGE SUMMARY')
print('======================================')
for col in ['Union_Coverage', 'Consensus_Coverage', 'Core_Coverage', 'ResourceSpecific_Coverage']:
    print(f'\n{col}')
    print('Mean   :', round(coverage_df[col].mean(), 4))
    print('Median :', round(coverage_df[col].median(), 4))
display(coverage_df[['Family_ID', 'Union_Genes', 'Union_Coverage', 'Consensus_Genes', 'Consensus_Coverage', 'Core_Genes', 'Core_Coverage', 'ResourceSpecific_Genes', 'ResourceSpecific_Coverage']])
OUT_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step3_family_coverage.csv'
coverage_df.to_csv(OUT_FILE, index=False)
print('\nSaved:')
print(OUT_FILE)
import pandas as pd
import numpy as np
rna_cols = [c for c in cptac_full.columns if c.endswith('_rna')]
cnv_cols = [c for c in cptac_full.columns if c.endswith('_cnv')]
X_rna = cptac_full[rna_cols].copy()
X_cnv = cptac_full[cnv_cols].copy()
X_rna_z = (X_rna - X_rna.mean(axis=0)) / X_rna.std(axis=0, ddof=0)
X_cnv_z = (X_cnv - X_cnv.mean(axis=0)) / X_cnv.std(axis=0, ddof=0)
X_rna_z = X_rna_z.replace([np.inf, -np.inf], np.nan)
X_cnv_z = X_cnv_z.replace([np.inf, -np.inf], np.nan)
activity = pd.DataFrame({'Cell_line': cptac_full['Cell_line'], 'Cancer_type': cptac_full['Cancer_type']})
feature_usage = []
for family_id, df in family_gene.groupby('Family_ID'):
    consensus_genes = sorted(set(df.loc[df['Resource_Support'] >= 2, 'Gene_Symbol']))
    rna_family_cols = [f'{g}_rna' for g in consensus_genes if f'{g}_rna' in X_rna_z.columns]
    cnv_family_cols = [f'{g}_cnv' for g in consensus_genes if f'{g}_cnv' in X_cnv_z.columns]
    if len(rna_family_cols) > 0:
        rna_score = X_rna_z[rna_family_cols].mean(axis=1, skipna=True)
    else:
        rna_score = pd.Series(np.nan, index=cptac_full.index)
    if len(cnv_family_cols) > 0:
        cnv_score = X_cnv_z[cnv_family_cols].mean(axis=1, skipna=True)
    else:
        cnv_score = pd.Series(np.nan, index=cptac_full.index)
    combined_score = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1, skipna=True)
    activity[family_id] = combined_score.values
    feature_usage.append({'Family_ID': family_id, 'Consensus_Genes_Total': len(consensus_genes), 'RNA_Genes_Used': len(rna_family_cols), 'CNV_Genes_Used': len(cnv_family_cols)})
feature_usage_df = pd.DataFrame(feature_usage)
family_cols = [c for c in activity.columns if c.startswith('BRCA_PRISM_F')]
print('======================================')
print('CPTAC CONSENSUS FAMILY ACTIVITY')
print('======================================')
print('\nPatients:')
print(len(activity))
print('\nFamilies:')
print(len(family_cols))
print('\nMissing family-score values:')
print(activity[family_cols].isna().sum().sum())
print('\nMean RNA genes used per family:')
print(round(feature_usage_df['RNA_Genes_Used'].mean(), 2))
print('\nMean CNV genes used per family:')
print(round(feature_usage_df['CNV_Genes_Used'].mean(), 2))
display(activity[['Cell_line', 'Cancer_type'] + family_cols[:5]].head())
display(feature_usage_df.head(10))
ACTIVITY_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step4_consensus_family_activity.csv'
FEATURE_FILE_OUT = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step4_features_used.csv'
activity.to_csv(ACTIVITY_FILE, index=False)
feature_usage_df.to_csv(FEATURE_FILE_OUT, index=False)
print('\nSaved:')
print(ACTIVITY_FILE)
print(FEATURE_FILE_OUT)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
ACTIVITY_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step4_consensus_family_activity.csv'
activity = pd.read_csv(ACTIVITY_FILE)
print('Activity shape:', activity.shape)
print('Columns:', activity.columns[:8].tolist())
display(activity.head())
import numpy as np
import pandas as pd
from scipy.stats import kruskal
family_cols = [c for c in activity.columns if c.startswith('BRCA_PRISM_F')]
subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
n = len(activity)
k = len(subtypes)
results = []
for family_id in family_cols:
    groups = [activity.loc[activity['Cancer_type'] == subtype, family_id].dropna().values for subtype in subtypes]
    H, p = kruskal(*groups)
    epsilon2 = (H - k + 1) / (n - k)
    epsilon2 = max(0, epsilon2)
    results.append({'Family_ID': family_id, 'Kruskal_H': H, 'P_Value': p, 'Epsilon_Squared': epsilon2})
cptac_assoc = pd.DataFrame(results)
pvals = cptac_assoc['P_Value'].to_numpy()
m = len(pvals)
order = np.argsort(pvals)
ranked_p = pvals[order]
adjusted = ranked_p * m / np.arange(1, m + 1)
adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
adjusted = np.minimum(adjusted, 1.0)
fdr = np.empty(m)
fdr[order] = adjusted
cptac_assoc['FDR'] = fdr
cptac_assoc['FDR_Significant'] = cptac_assoc['FDR'] < 0.05
cptac_assoc = cptac_assoc.sort_values('Epsilon_Squared', ascending=False).reset_index(drop=True)
cptac_assoc['Effect_Rank'] = np.arange(1, len(cptac_assoc) + 1)
print('======================================')
print('CPTAC PAM50 BIOLOGICAL REPLICATION')
print('======================================')
print('\nPatients:')
print(n)
print('\nFamilies tested:')
print(len(cptac_assoc))
print('\nFamilies significant at FDR < 0.05:')
print(cptac_assoc['FDR_Significant'].sum())
print('\nTop 10 families by effect size:')
display(cptac_assoc[['Family_ID', 'Epsilon_Squared', 'Kruskal_H', 'P_Value', 'FDR', 'Effect_Rank']].head(10))
OUT_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step5_PAM50_association.csv'
cptac_assoc.to_csv(OUT_FILE, index=False)
print('\nSaved:')
print(OUT_FILE)
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
TCGA_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_step9B_PAM50_family_association.csv'
CPTAC_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step5_PAM50_association.csv'
tcga_assoc = pd.read_csv(TCGA_FILE)
cptac_assoc = pd.read_csv(CPTAC_FILE)
print('TCGA file shape:', tcga_assoc.shape)
print('CPTAC file shape:', cptac_assoc.shape)
print('\nTCGA columns:')
print(tcga_assoc.columns.tolist())
print('\nCPTAC columns:')
print(cptac_assoc.columns.tolist())
possible_effect_cols = ['Epsilon_Squared', 'Epsilon2', 'Epsilon_Sq', 'epsilon_squared', 'epsilon2']
tcga_effect_col = None
for col in possible_effect_cols:
    if col in tcga_assoc.columns:
        tcga_effect_col = col
        break
if tcga_effect_col is None:
    raise ValueError('Could not identify TCGA epsilon-squared column.')
print('\nTCGA effect-size column:')
print(tcga_effect_col)

def normalize_family_id(x):
    x = str(x)
    if x.startswith('BRCA_PRISM_'):
        return x.replace('BRCA_PRISM_', '')
    return x
tcga_assoc['Family_ID_Normalized'] = tcga_assoc['Family_ID'].apply(normalize_family_id)
cptac_assoc['Family_ID_Normalized'] = cptac_assoc['Family_ID'].apply(normalize_family_id)
tcga_small = tcga_assoc[['Family_ID_Normalized', tcga_effect_col]].copy()
tcga_small = tcga_small.rename(columns={tcga_effect_col: 'TCGA_Epsilon_Squared'})
cptac_small = cptac_assoc[['Family_ID_Normalized', 'Epsilon_Squared', 'FDR']].copy()
cptac_small = cptac_small.rename(columns={'Epsilon_Squared': 'CPTAC_Epsilon_Squared', 'FDR': 'CPTAC_FDR'})
replication = tcga_small.merge(cptac_small, on='Family_ID_Normalized', how='inner')
print('\nFamilies matched:')
print(len(replication))
rho, p = spearmanr(replication['TCGA_Epsilon_Squared'], replication['CPTAC_Epsilon_Squared'])
print('\n======================================')
print('TCGA → CPTAC BIOLOGICAL REPLICATION')
print('======================================')
print('\nEffect-size Spearman correlation:')
print('rho =', round(rho, 4))
print('p   =', f'{p:.6e}')
replication['TCGA_Rank'] = replication['TCGA_Epsilon_Squared'].rank(ascending=False, method='min').astype(int)
replication['CPTAC_Rank'] = replication['CPTAC_Epsilon_Squared'].rank(ascending=False, method='min').astype(int)
tcga_top10 = set(replication.nsmallest(10, 'TCGA_Rank')['Family_ID_Normalized'])
cptac_top10 = set(replication.nsmallest(10, 'CPTAC_Rank')['Family_ID_Normalized'])
shared_top10 = tcga_top10 & cptac_top10
top10_jaccard = len(shared_top10) / len(tcga_top10 | cptac_top10)
print('\nTCGA Top-10 families:')
print(sorted(tcga_top10))
print('\nCPTAC Top-10 families:')
print(sorted(cptac_top10))
print('\nShared Top-10 families:')
print(sorted(shared_top10))
print('\nNumber shared:')
print(len(shared_top10), '/ 10')
print('\nTop-10 Jaccard:')
print(round(top10_jaccard, 4))
replication['Mean_Effect_Size'] = replication[['TCGA_Epsilon_Squared', 'CPTAC_Epsilon_Squared']].mean(axis=1)
replication = replication.sort_values('Mean_Effect_Size', ascending=False).reset_index(drop=True)
print('\nTop replicated families:')
display(replication[['Family_ID_Normalized', 'TCGA_Epsilon_Squared', 'CPTAC_Epsilon_Squared', 'TCGA_Rank', 'CPTAC_Rank', 'CPTAC_FDR']].head(15))
OUT_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step6_TCGA_CPTAC_effect_replication.csv'
replication.to_csv(OUT_FILE, index=False)
print('\nSaved:')
print(OUT_FILE)
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
family_cols = [c for c in activity.columns if c.startswith('BRCA_PRISM_F')]
cptac_medians = []
for family_id in family_cols:
    row = {'Family_ID': family_id.replace('BRCA_PRISM_', '')}
    for subtype in subtypes:
        row[f'CPTAC_{subtype}'] = activity.loc[activity['Cancer_type'] == subtype, family_id].median()
    cptac_medians.append(row)
cptac_medians = pd.DataFrame(cptac_medians)
tcga_pattern = tcga_assoc.copy()
tcga_pattern['Family_ID_Normalized'] = tcga_pattern['Family_ID'].astype(str).str.replace('BRCA_PRISM_', '', regex=False)
tcga_pattern = tcga_pattern.rename(columns={'Median_Basal': 'TCGA_Basal', 'Median_Her2': 'TCGA_Her2', 'Median_LumA': 'TCGA_LumA', 'Median_LumB': 'TCGA_LumB', 'Median_Normal-like': 'TCGA_Normal-like'})
pattern_df = tcga_pattern[['Family_ID_Normalized', 'TCGA_Basal', 'TCGA_Her2', 'TCGA_LumA', 'TCGA_LumB', 'TCGA_Normal-like']].merge(cptac_medians, left_on='Family_ID_Normalized', right_on='Family_ID', how='inner')
rows = []
for _, row in pattern_df.iterrows():
    tcga_values = np.array([row[f'TCGA_{s}'] for s in subtypes], dtype=float)
    cptac_values = np.array([row[f'CPTAC_{s}'] for s in subtypes], dtype=float)
    rho, rho_p = spearmanr(tcga_values, cptac_values)
    r, r_p = pearsonr(tcga_values, cptac_values)
    rows.append({'Family_ID': row['Family_ID_Normalized'], 'Subtype_Spearman': rho, 'Subtype_Spearman_P': rho_p, 'Subtype_Pearson': r, 'Subtype_Pearson_P': r_p})
pattern_concordance = pd.DataFrame(rows)
effect_small = replication[['Family_ID_Normalized', 'TCGA_Epsilon_Squared', 'CPTAC_Epsilon_Squared', 'TCGA_Rank', 'CPTAC_Rank']].copy()
effect_small = effect_small.rename(columns={'Family_ID_Normalized': 'Family_ID'})
pattern_concordance = pattern_concordance.merge(effect_small, on='Family_ID', how='left')
print('==========================================')
print('TCGA → CPTAC SUBTYPE-PATTERN CONCORDANCE')
print('==========================================')
print('\nMean subtype-pattern Spearman:')
print(round(pattern_concordance['Subtype_Spearman'].mean(), 4))
print('\nMedian subtype-pattern Spearman:')
print(round(pattern_concordance['Subtype_Spearman'].median(), 4))
print('\nFamilies with positive subtype-pattern correlation:')
print((pattern_concordance['Subtype_Spearman'] > 0).sum(), '/', len(pattern_concordance))
print('\nTop TCGA families and subtype-pattern replication:')
display(pattern_concordance.sort_values('TCGA_Rank')[['Family_ID', 'TCGA_Epsilon_Squared', 'CPTAC_Epsilon_Squared', 'TCGA_Rank', 'CPTAC_Rank', 'Subtype_Spearman', 'Subtype_Pearson']].head(15))
print('\nStrongest subtype-pattern agreement:')
display(pattern_concordance.sort_values('Subtype_Spearman', ascending=False)[['Family_ID', 'Subtype_Spearman', 'Subtype_Pearson', 'TCGA_Epsilon_Squared', 'CPTAC_Epsilon_Squared']].head(10))
OUT_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step7_subtype_pattern_concordance.csv'
pattern_concordance.to_csv(OUT_FILE, index=False)
print('\nSaved:')
print(OUT_FILE)
import pandas as pd
CPTAC_ASSOC_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step5_PAM50_association.csv'
cptac_assoc = pd.read_csv(CPTAC_ASSOC_FILE)
print(cptac_assoc[['Family_ID', 'Epsilon_Squared', 'FDR']].to_string(index=False))
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, classification_report
TCGA_DIR = PROJECT_ROOT / 'data_raw' / 'TCGA_BRCA'
CPTAC_DIR = PROJECT_ROOT / 'data_raw' / 'CPTAC'
TCGA_X_FILE = TCGA_DIR / 'tcga_brca_as_validation.csv'
TCGA_Y_FILE = TCGA_DIR / 'tcga_brca_mutation_cnv_rna_subtypes.csv'
CPTAC_X_FILE = CPTAC_DIR / 'cptac_as_validation.csv'
CPTAC_Y_FILE = CPTAC_DIR / 'cptac_brca_cnv_rna_subtypes_independent.csv'
FAMILY_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_step7_family_gene_support.csv'
tcga_x = pd.read_csv(TCGA_X_FILE)
tcga_y = pd.read_csv(TCGA_Y_FILE)
cptac_x = pd.read_csv(CPTAC_X_FILE)
cptac_y = pd.read_csv(CPTAC_Y_FILE)
family_gene = pd.read_csv(FAMILY_FILE)
tcga = tcga_x.merge(tcga_y[['Cell_line', 'Cancer_type']], on='Cell_line', how='inner')
cptac = cptac_x.merge(cptac_y[['Cell_line', 'Cancer_type']], on='Cell_line', how='inner')
print('TCGA samples :', len(tcga))
print('CPTAC samples:', len(cptac))
consensus = family_gene[family_gene['Resource_Support'] >= 2].copy()
families = sorted(consensus['Family_ID'].unique())
print('BRCA-PRISM families:', len(families))
all_consensus_genes = sorted(consensus['Gene_Symbol'].unique())
required_features = []
for gene in all_consensus_genes:
    rna = f'{gene}_rna'
    cnv = f'{gene}_cnv'
    if rna in tcga.columns and rna in cptac.columns:
        required_features.append(rna)
    if cnv in tcga.columns and cnv in cptac.columns:
        required_features.append(cnv)
required_features = sorted(set(required_features))
print('Common TCGA-CPTAC molecular features:', len(required_features))
tcga_features = tcga[required_features].copy()
cptac_features = cptac[required_features].copy()
tcga_medians = tcga_features.median(axis=0)
tcga_features = tcga_features.fillna(tcga_medians)
cptac_features = cptac_features.fillna(tcga_medians)
tcga_mean = tcga_features.mean(axis=0)
tcga_std = tcga_features.std(axis=0, ddof=0)
tcga_std = tcga_std.replace(0, 1.0)
tcga_z = (tcga_features - tcga_mean) / tcga_std
cptac_z = (cptac_features - tcga_mean) / tcga_std

def build_family_scores(z_df, consensus_df):
    scores = pd.DataFrame(index=z_df.index)
    usage = []
    for family_id, df in consensus_df.groupby('Family_ID'):
        genes = sorted(df['Gene_Symbol'].unique())
        rna_cols = [f'{g}_rna' for g in genes if f'{g}_rna' in z_df.columns]
        cnv_cols = [f'{g}_cnv' for g in genes if f'{g}_cnv' in z_df.columns]
        if rna_cols:
            rna_score = z_df[rna_cols].mean(axis=1)
        else:
            rna_score = pd.Series(np.nan, index=z_df.index)
        if cnv_cols:
            cnv_score = z_df[cnv_cols].mean(axis=1)
        else:
            cnv_score = pd.Series(np.nan, index=z_df.index)
        scores[family_id] = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1)
        usage.append({'Family_ID': family_id, 'RNA_Features': len(rna_cols), 'CNV_Features': len(cnv_cols)})
    return (scores, pd.DataFrame(usage))
X_tcga_family, usage = build_family_scores(tcga_z, consensus)
X_cptac_family, _ = build_family_scores(cptac_z, consensus)
X_cptac_family = X_cptac_family[X_tcga_family.columns]
print('\nTCGA family matrix:', X_tcga_family.shape)
print('CPTAC family matrix:', X_cptac_family.shape)
print('TCGA missing family scores:', X_tcga_family.isna().sum().sum())
print('CPTAC missing family scores:', X_cptac_family.isna().sum().sum())
y_tcga = tcga['Cancer_type'].values
y_cptac = cptac['Cancer_type'].values
print('\nTCGA class distribution:')
print(pd.Series(y_tcga).value_counts())
print('\nCPTAC class distribution:')
print(pd.Series(y_cptac).value_counts())
model = LogisticRegression(max_iter=5000, solver='lbfgs')
model.fit(X_tcga_family, y_tcga)
y_pred = model.predict(X_cptac_family)
accuracy = accuracy_score(y_cptac, y_pred)
balanced_acc = balanced_accuracy_score(y_cptac, y_pred)
macro_f1 = f1_score(y_cptac, y_pred, average='macro')
weighted_f1 = f1_score(y_cptac, y_pred, average='weighted')
mcc = matthews_corrcoef(y_cptac, y_pred)
print('\n========================================')
print('TCGA → CPTAC INDEPENDENT PAM50 TEST')
print('========================================')
print(f'Accuracy          : {accuracy:.4f}')
print(f'Balanced Accuracy : {balanced_acc:.4f}')
print(f'Macro-F1          : {macro_f1:.4f}')
print(f'Weighted-F1       : {weighted_f1:.4f}')
print(f'MCC               : {mcc:.4f}')
print('\nClassification report:\n')
print(classification_report(y_cptac, y_pred, digits=4, zero_division=0))
class_order = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
cm = confusion_matrix(y_cptac, y_pred, labels=class_order)
cm_df = pd.DataFrame(cm, index=class_order, columns=class_order)
print('\nConfusion Matrix:')
display(cm_df)
prediction_df = pd.DataFrame({'Cell_line': cptac['Cell_line'], 'True_PAM50': y_cptac, 'Predicted_PAM50': y_pred})
PRED_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step8_independent_predictions.csv'
prediction_df.to_csv(PRED_FILE, index=False)
metrics_df = pd.DataFrame({'Metric': ['Accuracy', 'Balanced Accuracy', 'Macro-F1', 'Weighted-F1', 'MCC'], 'Value': [accuracy, balanced_acc, macro_f1, weighted_f1, mcc]})
METRIC_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step8_independent_metrics.csv'
metrics_df.to_csv(METRIC_FILE, index=False)
print('\nSaved:')
print(PRED_FILE)
print(METRIC_FILE)
import numpy as np
import pandas as pd
print('======================================')
print('1. FAMILY-SCORE DISTRIBUTIONS')
print('======================================')
summary = pd.DataFrame({'TCGA_Mean': X_tcga_family.mean(), 'CPTAC_Mean': X_cptac_family.mean(), 'TCGA_SD': X_tcga_family.std(), 'CPTAC_SD': X_cptac_family.std()})
summary['Mean_Shift'] = summary['CPTAC_Mean'] - summary['TCGA_Mean']
summary['SD_Ratio'] = summary['CPTAC_SD'] / summary['TCGA_SD']
display(summary.sort_values('Mean_Shift', key=lambda x: x.abs(), ascending=False).head(15))
print('\n======================================')
print('2. OVERALL SCORE RANGE')
print('======================================')
print('TCGA min/max:', X_tcga_family.min().min(), X_tcga_family.max().max())
print('CPTAC min/max:', X_cptac_family.min().min(), X_cptac_family.max().max())
print('\nMean absolute TCGA family score:', X_tcga_family.abs().mean().mean())
print('Mean absolute CPTAC family score:', X_cptac_family.abs().mean().mean())
print('\n======================================')
print('3. PREDICTED CLASS COUNTS')
print('======================================')
print(pd.Series(y_pred).value_counts())
print('\n======================================')
print('4. PREDICTED PROBABILITIES')
print('======================================')
proba = model.predict_proba(X_cptac_family)
proba_df = pd.DataFrame(proba, columns=model.classes_)
display(proba_df.describe().T)
print('\n======================================')
print('5. FIRST 10 CPTAC PREDICTIONS')
print('======================================')
check_df = pd.DataFrame({'True': y_cptac, 'Predicted': y_pred})
for cls in model.classes_:
    check_df[f'P_{cls}'] = proba_df[cls].values
display(check_df.head(10))
print('\n======================================')
print('6. MODEL INTERCEPTS')
print('======================================')
print(pd.Series(model.intercept_, index=model.classes_))
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, classification_report
tcga_raw = tcga[required_features].copy()
cptac_raw = cptac[required_features].copy()
tcga_medians_local = tcga_raw.median(axis=0)
cptac_medians_local = cptac_raw.median(axis=0)
tcga_imp = tcga_raw.fillna(tcga_medians_local)
cptac_imp = cptac_raw.fillna(cptac_medians_local)
tcga_mean_local = tcga_imp.mean(axis=0)
tcga_std_local = tcga_imp.std(axis=0, ddof=0).replace(0, 1.0)
cptac_mean_local = cptac_imp.mean(axis=0)
cptac_std_local = cptac_imp.std(axis=0, ddof=0).replace(0, 1.0)
tcga_z_local = (tcga_imp - tcga_mean_local) / tcga_std_local
cptac_z_local = (cptac_imp - cptac_mean_local) / cptac_std_local

def build_consensus_family_scores(z_df, consensus_df):
    scores = pd.DataFrame(index=z_df.index)
    for family_id, df in consensus_df.groupby('Family_ID'):
        genes = sorted(df['Gene_Symbol'].unique())
        rna_cols = [f'{g}_rna' for g in genes if f'{g}_rna' in z_df.columns]
        cnv_cols = [f'{g}_cnv' for g in genes if f'{g}_cnv' in z_df.columns]
        if len(rna_cols) > 0:
            rna_score = z_df[rna_cols].mean(axis=1)
        else:
            rna_score = pd.Series(np.nan, index=z_df.index)
        if len(cnv_cols) > 0:
            cnv_score = z_df[cnv_cols].mean(axis=1)
        else:
            cnv_score = pd.Series(np.nan, index=z_df.index)
        scores[family_id] = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1)
    return scores
X_tcga_local = build_consensus_family_scores(tcga_z_local, consensus)
X_cptac_local = build_consensus_family_scores(cptac_z_local, consensus)
X_cptac_local = X_cptac_local[X_tcga_local.columns]
print('======================================')
print('NORMALIZED FAMILY-SCORE CHECK')
print('======================================')
print('\nTCGA family matrix:')
print(X_tcga_local.shape)
print('CPTAC family matrix:')
print(X_cptac_local.shape)
print('\nMean absolute TCGA family score:')
print(X_tcga_local.abs().mean().mean())
print('\nMean absolute CPTAC family score:')
print(X_cptac_local.abs().mean().mean())
print('\nTCGA range:')
print(X_tcga_local.min().min(), X_tcga_local.max().max())
print('\nCPTAC range:')
print(X_cptac_local.min().min(), X_cptac_local.max().max())
model_local = LogisticRegression(max_iter=5000, solver='lbfgs')
model_local.fit(X_tcga_local, y_tcga)
y_pred_local = model_local.predict(X_cptac_local)
y_prob_local = model_local.predict_proba(X_cptac_local)
accuracy = accuracy_score(y_cptac, y_pred_local)
balanced_acc = balanced_accuracy_score(y_cptac, y_pred_local)
macro_f1 = f1_score(y_cptac, y_pred_local, average='macro')
weighted_f1 = f1_score(y_cptac, y_pred_local, average='weighted')
mcc = matthews_corrcoef(y_cptac, y_pred_local)
print('\n======================================')
print('TCGA -> CPTAC SENSITIVITY ANALYSIS')
print('Cohort-wise unsupervised normalization')
print('======================================')
print(f'\nAccuracy          : {accuracy:.4f}')
print(f'Balanced Accuracy : {balanced_acc:.4f}')
print(f'Macro-F1          : {macro_f1:.4f}')
print(f'Weighted-F1       : {weighted_f1:.4f}')
print(f'MCC               : {mcc:.4f}')
print('\nPredicted class counts:')
print(pd.Series(y_pred_local).value_counts())
print('\nClassification report:\n')
print(classification_report(y_cptac, y_pred_local, digits=4, zero_division=0))
class_order = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
cm = confusion_matrix(y_cptac, y_pred_local, labels=class_order)
cm_df = pd.DataFrame(cm, index=class_order, columns=class_order)
print('\nConfusion Matrix:')
display(cm_df)
prediction_df = pd.DataFrame({'Cell_line': cptac['Cell_line'], 'True_PAM50': y_cptac, 'Predicted_PAM50': y_pred_local})
PRED_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step8C_cohort_normalized_predictions.csv'
prediction_df.to_csv(PRED_FILE, index=False)
metrics_df = pd.DataFrame({'Metric': ['Accuracy', 'Balanced Accuracy', 'Macro-F1', 'Weighted-F1', 'MCC'], 'Value': [accuracy, balanced_acc, macro_f1, weighted_f1, mcc]})
METRIC_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step8C_cohort_normalized_metrics.csv'
metrics_df.to_csv(METRIC_FILE, index=False)
print('\nSaved:')
print(PRED_FILE)
print(METRIC_FILE)
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef, confusion_matrix

def build_ablation_family_scores(z_df, family_gene_df, mode):
    scores = pd.DataFrame(index=z_df.index)
    for family_id, df in family_gene_df.groupby('Family_ID'):
        if mode == 'All Family Genes':
            genes = sorted(df['Gene_Symbol'].unique())
        elif mode == 'Consensus Genes':
            genes = sorted(df.loc[df['Resource_Support'] >= 2, 'Gene_Symbol'].unique())
        elif mode == 'Resource-Specific Genes':
            genes = sorted(df.loc[df['Resource_Support'] == 1, 'Gene_Symbol'].unique())
        else:
            raise ValueError(f'Unknown mode: {mode}')
        rna_cols = [f'{g}_rna' for g in genes if f'{g}_rna' in z_df.columns]
        cnv_cols = [f'{g}_cnv' for g in genes if f'{g}_cnv' in z_df.columns]
        if len(rna_cols) == 0 and len(cnv_cols) == 0:
            continue
        if len(rna_cols) > 0:
            rna_score = z_df[rna_cols].mean(axis=1)
        else:
            rna_score = pd.Series(np.nan, index=z_df.index)
        if len(cnv_cols) > 0:
            cnv_score = z_df[cnv_cols].mean(axis=1)
        else:
            cnv_score = pd.Series(np.nan, index=z_df.index)
        scores[family_id] = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1, skipna=True)
    return scores
print('TCGA normalized molecular matrix:', tcga_z_local.shape)
print('CPTAC normalized molecular matrix:', cptac_z_local.shape)
modes = ['All Family Genes', 'Consensus Genes', 'Resource-Specific Genes']
results = []
prediction_tables = {}
confusion_tables = {}
for mode in modes:
    print('\n======================================')
    print(mode)
    print('======================================')
    X_tcga_mode = build_ablation_family_scores(tcga_z_local, family_gene, mode)
    X_cptac_mode = build_ablation_family_scores(cptac_z_local, family_gene, mode)
    common_families = sorted(set(X_tcga_mode.columns) & set(X_cptac_mode.columns))
    X_tcga_mode = X_tcga_mode[common_families]
    X_cptac_mode = X_cptac_mode[common_families]
    print('Families used:', len(common_families))
    print('TCGA missing:', X_tcga_mode.isna().sum().sum())
    print('CPTAC missing:', X_cptac_mode.isna().sum().sum())
    X_tcga_mode = X_tcga_mode.fillna(0)
    X_cptac_mode = X_cptac_mode.fillna(0)
    clf = LogisticRegression(max_iter=5000, solver='lbfgs')
    clf.fit(X_tcga_mode, y_tcga)
    pred = clf.predict(X_cptac_mode)
    acc = accuracy_score(y_cptac, pred)
    bal = balanced_accuracy_score(y_cptac, pred)
    macro = f1_score(y_cptac, pred, average='macro')
    weighted = f1_score(y_cptac, pred, average='weighted')
    mcc = matthews_corrcoef(y_cptac, pred)
    results.append({'Representation': mode, 'Families_Used': len(common_families), 'Accuracy': acc, 'Balanced_Accuracy': bal, 'Macro_F1': macro, 'Weighted_F1': weighted, 'MCC': mcc})
    print(f'Accuracy          : {acc:.4f}')
    print(f'Balanced Accuracy : {bal:.4f}')
    print(f'Macro-F1          : {macro:.4f}')
    print(f'Weighted-F1       : {weighted:.4f}')
    print(f'MCC               : {mcc:.4f}')
    print('\nPredicted classes:')
    print(pd.Series(pred).value_counts())
    prediction_tables[mode] = pd.DataFrame({'Cell_line': cptac['Cell_line'], 'True_PAM50': y_cptac, 'Predicted_PAM50': pred})
    class_order = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
    cm = confusion_matrix(y_cptac, pred, labels=class_order)
    confusion_tables[mode] = pd.DataFrame(cm, index=class_order, columns=class_order)
ablation_external = pd.DataFrame(results)
ablation_external = ablation_external.sort_values('Macro_F1', ascending=False).reset_index(drop=True)
print('\n======================================')
print('EXTERNAL ABLATION SUMMARY')
print('======================================')
display(ablation_external)
for mode in modes:
    print('\n', mode)
    display(confusion_tables[mode])
SUMMARY_FILE = PROJECT_ROOT / 'results' / 'BRCA_PRISM_CPTAC_step9_external_ablation_summary.csv'
ablation_external.to_csv(SUMMARY_FILE, index=False)
for mode, df in prediction_tables.items():
    safe_name = mode.lower().replace(' ', '_').replace('-', '_')
    out = PROJECT_ROOT / 'results' / f'BRCA_PRISM_CPTAC_step9_{safe_name}_predictions.csv'
    df.to_csv(out, index=False)
print('\nSaved:')
print(SUMMARY_FILE)
import os
import pandas as pd
for f in sorted(os.listdir('results')):
    if f.endswith('.csv') or f.endswith('.xls'):
        try:
            cols = [c.lower() for c in pd.read_csv(f'results/{f}', nrows=0).columns]
            has_family = any(('family' in c and 'id' in c for c in cols))
            has_gene = any((c == 'gene_symbol' or c == 'gene' for c in cols))
            if has_family and has_gene:
                print(f)
        except Exception:
            pass
import os
import pandas as pd
for f in sorted(os.listdir('results')):
    if f.endswith('.csv') or f.endswith('.xls'):
        try:
            cols = [c.lower() for c in pd.read_csv(f'results/{f}', nrows=0).columns]
            has_family = any(('family' in c for c in cols))
            has_gene = any(('gene' in c for c in cols))
            if has_family and has_gene:
                df = pd.read_csv(f'results/{f}')
                fam_col = [c for c in df.columns if 'family' in c.lower()][0]
                if 'F0025' in str(df[fam_col].unique()[:5]) or df[fam_col].astype(str).str.contains('F0025').any():
                    print(f, '-> HAS F0025 rows')
        except Exception:
            pass
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import confusion_matrix
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
REPLICATION_FILE = RESULTS_DIR / 'BRCA_PRISM_CPTAC_step6_TCGA_CPTAC_effect_replication.csv'
PATTERN_FILE = RESULTS_DIR / 'BRCA_PRISM_CPTAC_step7_subtype_pattern_concordance.csv'
ABLATION_FILE = RESULTS_DIR / 'BRCA_PRISM_CPTAC_step9_external_ablation_summary.csv'
CONSENSUS_PRED_FILE = RESULTS_DIR / 'BRCA_PRISM_CPTAC_step9_consensus_genes_predictions.csv'
replication = pd.read_csv(REPLICATION_FILE)
pattern = pd.read_csv(PATTERN_FILE)
ablation = pd.read_csv(ABLATION_FILE)
predictions = pd.read_csv(CONSENSUS_PRED_FILE)
print('Replication:', replication.shape)
print('Pattern:', pattern.shape)
print('Ablation:', ablation.shape)
print('Predictions:', predictions.shape)
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
ax = axes[0, 0]
x = replication['TCGA_Epsilon_Squared']
y = replication['CPTAC_Epsilon_Squared']
ax.scatter(x, y, s=55, alpha=0.75, edgecolor='black', linewidth=0.5)
max_val = max(x.max(), y.max()) + 0.03
ax.plot([0, max_val], [0, max_val], linestyle='--', linewidth=1.2, color='gray')
rho_effect, p_effect = spearmanr(x, y)
highlights = {'F0006': 'F0006\nDNA replication', 'F0015': 'F0015\nCell cycle', 'F0010': 'F0010\nMismatch repair', 'F0023': 'F0023'}
for family_id, label in highlights.items():
    row = replication.loc[replication['Family_ID_Normalized'] == family_id]
    if len(row) == 1:
        xx = row['TCGA_Epsilon_Squared'].iloc[0]
        yy = row['CPTAC_Epsilon_Squared'].iloc[0]
        ax.scatter(xx, yy, s=110, edgecolor='black', linewidth=1.2)
        ax.annotate(label, (xx, yy), xytext=(7, 7), textcoords='offset points', fontsize=8.5, fontweight='bold')
ax.text(0.04, 0.95, f'All 31 families\nSpearman ρ = {rho_effect:.3f}\np = {p_effect:.3f}', transform=ax.transAxes, ha='left', va='top', fontsize=10, bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.7'))
ax.set_xlabel('TCGA PAM50 Effect Size (ε²)', fontsize=11)
ax.set_ylabel('CPTAC PAM50 Effect Size (ε²)', fontsize=11)
ax.set_xlim(-0.01, max_val)
ax.set_ylim(-0.01, max_val)
ax.grid(alpha=0.25, linestyle='--')
ax.set_title('A. Cross-Cohort Pathway Effect Sizes', fontsize=13, fontweight='bold')
ax = axes[0, 1]
pattern_plot = pattern[['Family_ID', 'Subtype_Spearman']].sort_values('Subtype_Spearman', ascending=True).reset_index(drop=True)
y_pos = np.arange(len(pattern_plot))
ax.barh(y_pos, pattern_plot['Subtype_Spearman'], height=0.72)
ax.axvline(0, linewidth=1, color='black')
mean_rho = pattern_plot['Subtype_Spearman'].mean()
median_rho = pattern_plot['Subtype_Spearman'].median()
positive_n = (pattern_plot['Subtype_Spearman'] > 0).sum()
ax.axvline(mean_rho, linestyle='--', linewidth=1.5, color='gray')
ax.set_yticks(y_pos)
ax.set_yticklabels(pattern_plot['Family_ID'], fontsize=7)
ax.set_xlim(-1.05, 1.05)
ax.set_xlabel('TCGA–CPTAC Subtype-Pattern Spearman ρ', fontsize=11)
ax.text(0.03, 0.97, f'Mean ρ = {mean_rho:.3f}\nMedian ρ = {median_rho:.2f}\nPositive = {positive_n}/31', transform=ax.transAxes, ha='left', va='top', fontsize=10, bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.7'))
ax.grid(axis='x', alpha=0.25, linestyle='--')
ax.set_title('B. PAM50 Subtype-Pattern Concordance', fontsize=13, fontweight='bold')
ax = axes[1, 0]
class_order = ['Basal', 'Her2', 'LumA', 'LumB', 'Normal-like']
cm = confusion_matrix(predictions['True_PAM50'], predictions['Predicted_PAM50'], labels=class_order)
row_sum = cm.sum(axis=1, keepdims=True)
cm_pct = np.divide(cm, row_sum, out=np.zeros_like(cm, dtype=float), where=row_sum != 0) * 100
im = ax.imshow(cm_pct, vmin=0, vmax=100, cmap='Blues')
ax.set_xticks(np.arange(len(class_order)))
ax.set_yticks(np.arange(len(class_order)))
ax.set_xticklabels(class_order, rotation=35, ha='right', fontsize=9)
ax.set_yticklabels(class_order, fontsize=9)
ax.set_xlabel('Predicted PAM50', fontsize=11)
ax.set_ylabel('True PAM50', fontsize=11)
for i in range(len(class_order)):
    for j in range(len(class_order)):
        text_color = 'white' if cm_pct[i, j] > 55 else 'black'
        ax.text(j, i, f'{cm[i, j]}\n({cm_pct[i, j]:.0f}%)', ha='center', va='center', fontsize=9, color=text_color)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Row Percentage (%)', fontsize=9)
ax.set_title('C. TCGA-Trained Consensus Model → CPTAC', fontsize=13, fontweight='bold')
ax.text(0.02, -0.27, 'Logistic Regression | 31 consensus families\nAccuracy = 0.656 | Balanced Accuracy = 0.442 | Macro-F1 = 0.447 | MCC = 0.481', transform=ax.transAxes, fontsize=9.5, ha='left', va='top')
ax = axes[1, 1]
representation_order = ['All Family Genes', 'Consensus Genes', 'Resource-Specific Genes']
ablation_plot = ablation.set_index('Representation').loc[representation_order].reset_index()
metrics = ['Accuracy', 'Balanced_Accuracy', 'Macro_F1', 'MCC']
metric_labels = ['Accuracy', 'Balanced\nAccuracy', 'Macro-F1', 'MCC']
x_pos = np.arange(len(metrics))
width = 0.24
for i, row in ablation_plot.iterrows():
    values = [row[m] for m in metrics]
    positions = x_pos + (i - 1) * width
    bars = ax.bar(positions, values, width=width, label=row['Representation'])
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.012, f'{value:.3f}', ha='center', va='bottom', fontsize=7.5, rotation=90)
ax.set_xticks(x_pos)
ax.set_xticklabels(metric_labels, fontsize=9)
ax.set_ylabel('Performance', fontsize=11)
ax.set_ylim(0, 0.78)
ax.grid(axis='y', alpha=0.25, linestyle='--')
ax.legend(loc='upper right', frameon=False, fontsize=8)
ax.text(0.5, -0.16, 'Families used: All = 31 | Consensus = 31 | Resource-specific = 27', transform=ax.transAxes, ha='center', va='top', fontsize=8.5)
ax.set_title('D. External Gene-Layer Ablation', fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0.04, 1, 1], h_pad=4.0, w_pad=2.5)
PNG_FILE = FIGURES_DIR / 'BRCA_PRISM_CPTAC_external_validation_composite.png'
PDF_FILE = FIGURES_DIR / 'BRCA_PRISM_CPTAC_external_validation_composite.pdf'
plt.savefig(PNG_FILE, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.savefig(PDF_FILE, bbox_inches='tight', pad_inches=0.2)
plt.show()
print('\nSaved:')
print(PNG_FILE)
print(PDF_FILE)
