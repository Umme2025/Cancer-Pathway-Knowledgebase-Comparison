from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
GENE_SUPPORT_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
tcga = pd.read_csv(TCGA_FILE)
gene_support = pd.read_csv(GENE_SUPPORT_FILE)
gene_support['Gene_Symbol'] = gene_support['Gene_Symbol'].astype(str).str.upper().str.strip()
rna_cols = [c for c in tcga.columns if c.lower().endswith('_rna')]
cnv_cols = [c for c in tcga.columns if c.lower().endswith('_cnv')]
feature_cols = rna_cols + cnv_cols
Xraw = tcga[feature_cols].apply(pd.to_numeric, errors='coerce')
means = Xraw.mean(axis=0)
stds = Xraw.std(axis=0, ddof=0).replace(0, np.nan)
Z = (Xraw - means) / stds

def build_family_scores(mode):
    score_df = tcga[['Cell_line', 'PAM50_Subtype']].copy()
    for family_id, g in gene_support.groupby('Family_ID'):
        if mode == 'Consensus':
            genes = set(g.loc[g['Resource_Support'] >= 2, 'Gene_Symbol'])
        elif mode == 'Resource_Specific':
            genes = set(g.loc[g['Resource_Support'] == 1, 'Gene_Symbol'])
        elif mode == 'All':
            genes = set(g['Gene_Symbol'])
        fam_rna = [c for c in rna_cols if c[:-4].upper() in genes]
        fam_cnv = [c for c in cnv_cols if c[:-4].upper() in genes]
        if len(fam_rna) == 0 and len(fam_cnv) == 0:
            score_df[family_id] = np.nan
            continue
        rna_score = Z[fam_rna].mean(axis=1) if fam_rna else pd.Series(np.nan, index=Z.index)
        cnv_score = Z[fam_cnv].mean(axis=1) if fam_cnv else pd.Series(np.nan, index=Z.index)
        combined = pd.concat([rna_score, cnv_score], axis=1).mean(axis=1)
        score_df[family_id] = combined
    return score_df
all_scores = build_family_scores('All')
consensus_scores = build_family_scores('Consensus')
specific_scores = build_family_scores('Resource_Specific')
family_cols = [c for c in all_scores.columns if c.startswith('BRCA_PRISM_F')]
common_families = []
for fam in family_cols:
    if all_scores[fam].notna().any() and consensus_scores[fam].notna().any() and specific_scores[fam].notna().any():
        common_families.append(fam)
print('Families available in all ablations:', len(common_families))
datasets = {'All_Family_Genes': all_scores[common_families].copy(), 'Consensus_Genes': consensus_scores[common_families].copy(), 'Resource_Specific_Genes': specific_scores[common_families].copy()}
y = tcga['PAM50_Subtype'].copy()
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = []
for representation, X in datasets.items():
    X = X.fillna(X.mean())
    print('\nRunning:', representation, 'features:', X.shape[1])
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        results.append({'Representation': representation, 'Fold': fold, 'Number_of_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Number_of_Families=('Number_of_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 16 CONSENSUS VS RESOURCE-SPECIFIC ABLATION')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16_ablation_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16_ablation_summary.csv', index=False)
print('\nSTEP 16 COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
GENE_SUPPORT_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
tcga = pd.read_csv(TCGA_FILE)
gene_support = pd.read_csv(GENE_SUPPORT_FILE)
gene_support['Gene_Symbol'] = gene_support['Gene_Symbol'].astype(str).str.upper().str.strip()
rna_cols = [c for c in tcga.columns if c.lower().endswith('_rna')]
cnv_cols = [c for c in tcga.columns if c.lower().endswith('_cnv')]
feature_cols = rna_cols + cnv_cols
Xraw = tcga[feature_cols].apply(pd.to_numeric, errors='coerce')
y = tcga['PAM50_Subtype'].copy()
family_gene_sets = {}
for family_id, g in gene_support.groupby('Family_ID'):
    all_genes = set(g['Gene_Symbol'])
    consensus_genes = set(g.loc[g['Resource_Support'] >= 2, 'Gene_Symbol'])
    specific_genes = set(g.loc[g['Resource_Support'] == 1, 'Gene_Symbol'])
    family_gene_sets[family_id] = {'All_Family_Genes': all_genes, 'Consensus_Genes': consensus_genes, 'Resource_Specific_Genes': specific_genes}

def measurable_gene_count(genes):
    rna = sum((1 for g in genes if f'{g}_rna' in Xraw.columns))
    cnv = sum((1 for g in genes if f'{g}_cnv' in Xraw.columns))
    return rna + cnv
common_families = []
for family_id, modes in family_gene_sets.items():
    counts = [measurable_gene_count(modes['All_Family_Genes']), measurable_gene_count(modes['Consensus_Genes']), measurable_gene_count(modes['Resource_Specific_Genes'])]
    if all((c > 0 for c in counts)):
        common_families.append(family_id)
common_families = sorted(common_families)
print('Families available in all ablations:', len(common_families))
gene_count_rows = []
for family_id in common_families:
    for mode in ['All_Family_Genes', 'Consensus_Genes', 'Resource_Specific_Genes']:
        genes = family_gene_sets[family_id][mode]
        rna_n = sum((1 for g in genes if f'{g}_rna' in Xraw.columns))
        cnv_n = sum((1 for g in genes if f'{g}_cnv' in Xraw.columns))
        gene_count_rows.append({'Family_ID': family_id, 'Representation': mode, 'Unique_Genes': len(genes), 'RNA_Features': rna_n, 'CNV_Features': cnv_n})
gene_counts = pd.DataFrame(gene_count_rows)
print('\nMean gene counts per family:')
display(gene_counts.groupby('Representation').agg(Mean_Unique_Genes=('Unique_Genes', 'mean'), Median_Unique_Genes=('Unique_Genes', 'median'), Mean_RNA_Features=('RNA_Features', 'mean'), Mean_CNV_Features=('CNV_Features', 'mean')).round(2))
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
splits = list(cv.split(Xraw, y))
results = []
for mode in ['All_Family_Genes', 'Consensus_Genes', 'Resource_Specific_Genes']:
    print('\nRunning:', mode)
    for fold, (train_idx, test_idx) in enumerate(splits, start=1):
        X_train_raw = Xraw.iloc[train_idx].copy()
        X_test_raw = Xraw.iloc[test_idx].copy()
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        train_means = X_train_raw.mean(axis=0)
        train_stds = X_train_raw.std(axis=0, ddof=0).replace(0, np.nan)
        Z_train = (X_train_raw - train_means) / train_stds
        Z_test = (X_test_raw - train_means) / train_stds
        train_family_scores = pd.DataFrame(index=X_train_raw.index)
        test_family_scores = pd.DataFrame(index=X_test_raw.index)
        for family_id in common_families:
            genes = family_gene_sets[family_id][mode]
            fam_rna = [f'{g}_rna' for g in genes if f'{g}_rna' in Xraw.columns]
            fam_cnv = [f'{g}_cnv' for g in genes if f'{g}_cnv' in Xraw.columns]
            if fam_rna:
                train_rna = Z_train[fam_rna].mean(axis=1)
                test_rna = Z_test[fam_rna].mean(axis=1)
            else:
                train_rna = pd.Series(np.nan, index=X_train_raw.index)
                test_rna = pd.Series(np.nan, index=X_test_raw.index)
            if fam_cnv:
                train_cnv = Z_train[fam_cnv].mean(axis=1)
                test_cnv = Z_test[fam_cnv].mean(axis=1)
            else:
                train_cnv = pd.Series(np.nan, index=X_train_raw.index)
                test_cnv = pd.Series(np.nan, index=X_test_raw.index)
            train_combined = pd.concat([train_rna, train_cnv], axis=1).mean(axis=1)
            test_combined = pd.concat([test_rna, test_cnv], axis=1).mean(axis=1)
            train_family_scores[family_id] = train_combined
            test_family_scores[family_id] = test_combined
        family_means = train_family_scores.mean()
        train_family_scores = train_family_scores.fillna(family_means)
        test_family_scores = test_family_scores.fillna(family_means)
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(train_family_scores, y_train)
        pred = model.predict(test_family_scores)
        results.append({'Representation': mode, 'Fold': fold, 'Number_of_Families': len(common_families), 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Number_of_Families=('Number_of_Families', 'first'), Mean_Accuracy=('Accuracy', 'mean'), SD_Accuracy=('Accuracy', 'std'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), SD_Balanced_Accuracy=('Balanced_Accuracy', 'std'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 16B LEAKAGE-CORRECTED ABLATION')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16B_corrected_ablation_cv_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16B_corrected_ablation_summary.csv', index=False)
gene_counts.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16B_ablation_gene_counts.csv', index=False)
print('\nSTEP 16B COMPLETE')
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = RESULTS_DIR / 'BRCA_PRISM_step9A_TCGA_with_PAM50_labels.csv'
GENE_SUPPORT_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
tcga = pd.read_csv(TCGA_FILE)
gene_support = pd.read_csv(GENE_SUPPORT_FILE)
gene_support['Gene_Symbol'] = gene_support['Gene_Symbol'].astype(str).str.upper().str.strip()
rna_cols = [c for c in tcga.columns if c.lower().endswith('_rna')]
cnv_cols = [c for c in tcga.columns if c.lower().endswith('_cnv')]
feature_cols = rna_cols + cnv_cols
Xraw = tcga[feature_cols].apply(pd.to_numeric, errors='coerce')
y = tcga['PAM50_Subtype'].copy()
family_gene_sets = {}
for family_id, g in gene_support.groupby('Family_ID'):
    family_gene_sets[family_id] = {'All_Family_Genes': set(g['Gene_Symbol']), 'Consensus_Genes': set(g.loc[g['Resource_Support'] >= 2, 'Gene_Symbol']), 'Resource_Specific_Genes': set(g.loc[g['Resource_Support'] == 1, 'Gene_Symbol'])}

def measurable(genes):
    return any((f'{g}_rna' in Xraw.columns or f'{g}_cnv' in Xraw.columns for g in genes))
common_families = []
for family_id, modes in family_gene_sets.items():
    if all((measurable(modes[m]) for m in ['All_Family_Genes', 'Consensus_Genes', 'Resource_Specific_Genes'])):
        common_families.append(family_id)
common_families = sorted(common_families)
print('Families used:', len(common_families))
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
splits = list(cv.split(Xraw, y))
print('Total CV evaluations:', len(splits))
results = []
for mode in ['All_Family_Genes', 'Consensus_Genes', 'Resource_Specific_Genes']:
    print('\nRunning:', mode)
    for split_no, (train_idx, test_idx) in enumerate(splits, start=1):
        X_train_raw = Xraw.iloc[train_idx]
        X_test_raw = Xraw.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]
        mu = X_train_raw.mean()
        sd = X_train_raw.std(ddof=0).replace(0, np.nan)
        Z_train = (X_train_raw - mu) / sd
        Z_test = (X_test_raw - mu) / sd
        train_scores = pd.DataFrame(index=train_idx)
        test_scores = pd.DataFrame(index=test_idx)
        for family_id in common_families:
            genes = family_gene_sets[family_id][mode]
            fam_rna = [f'{g}_rna' for g in genes if f'{g}_rna' in Xraw.columns]
            fam_cnv = [f'{g}_cnv' for g in genes if f'{g}_cnv' in Xraw.columns]
            tr_parts = []
            te_parts = []
            if fam_rna:
                tr_parts.append(Z_train[fam_rna].mean(axis=1))
                te_parts.append(Z_test[fam_rna].mean(axis=1))
            if fam_cnv:
                tr_parts.append(Z_train[fam_cnv].mean(axis=1))
                te_parts.append(Z_test[fam_cnv].mean(axis=1))
            train_scores[family_id] = pd.concat(tr_parts, axis=1).mean(axis=1)
            test_scores[family_id] = pd.concat(te_parts, axis=1).mean(axis=1)
        train_means = train_scores.mean()
        train_scores = train_scores.fillna(train_means)
        test_scores = test_scores.fillna(train_means)
        model = LogisticRegression(max_iter=3000, class_weight='balanced', solver='lbfgs')
        model.fit(train_scores, y_train)
        pred = model.predict(test_scores)
        repeat_id = (split_no - 1) // 5 + 1
        fold_id = (split_no - 1) % 5 + 1
        results.append({'Representation': mode, 'Repeat': repeat_id, 'Fold': fold_id, 'Accuracy': accuracy_score(y_test, pred), 'Balanced_Accuracy': balanced_accuracy_score(y_test, pred), 'Macro_F1': f1_score(y_test, pred, average='macro')})
results = pd.DataFrame(results)
summary = results.groupby('Representation').agg(Mean_Accuracy=('Accuracy', 'mean'), SD_Accuracy=('Accuracy', 'std'), Mean_Balanced_Accuracy=('Balanced_Accuracy', 'mean'), SD_Balanced_Accuracy=('Balanced_Accuracy', 'std'), Mean_Macro_F1=('Macro_F1', 'mean'), SD_Macro_F1=('Macro_F1', 'std')).reset_index().sort_values('Mean_Macro_F1', ascending=False)
print('\n' + '=' * 100)
print('STEP 16C REPEATED-CV ABLATION')
print('=' * 100)
display(summary.round(4))
results.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16C_repeatedCV_ablation_results.csv', index=False)
summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step16C_repeatedCV_ablation_summary.csv', index=False)
print('\nSTEP 16C COMPLETE')
from pathlib import Path
import pandas as pd
from scipy.stats import wilcoxon
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FILE = RESULTS_DIR / 'BRCA_PRISM_step16C_repeatedCV_ablation_results.csv'
df = pd.read_csv(FILE)
pivot = df.pivot_table(index=['Repeat', 'Fold'], columns='Representation', values='Macro_F1').reset_index()
print('Paired CV splits:', len(pivot))
comparisons = [('Resource_Specific_Genes', 'Consensus_Genes'), ('Resource_Specific_Genes', 'All_Family_Genes'), ('All_Family_Genes', 'Consensus_Genes')]
rows = []
for a, b in comparisons:
    diff = pivot[a] - pivot[b]
    stat, p = wilcoxon(pivot[a], pivot[b], alternative='two-sided')
    rows.append({'Comparison': f'{a} vs {b}', 'Mean_Macro_F1_A': pivot[a].mean(), 'Mean_Macro_F1_B': pivot[b].mean(), 'Mean_Paired_Difference': diff.mean(), 'Median_Paired_Difference': diff.median(), 'Wilcoxon_Statistic': stat, 'P_Value': p})
results = pd.DataFrame(rows)
results = results.sort_values('P_Value').reset_index(drop=True)
m = len(results)
results['FDR'] = results['P_Value'] * m / (results.index + 1)
results['FDR'] = results['FDR'][::-1].cummin()[::-1].clip(upper=1)
print('\n' + '=' * 100)
print('STEP 16D PAIRED STATISTICAL TEST')
print('=' * 100)
display(results.round(6))
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step16D_paired_ablation_statistics.csv'
results.to_csv(OUT_FILE, index=False)
print('\nSTEP 16D COMPLETE')
print('Saved:', OUT_FILE)
