from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import matplotlib.pyplot as plt
FAMILIES = {'F0006': 'DNA replication', 'F0025': 'Endometrial cancer'}
BRCA_EFFECT_SIZE = {'F0006': 0.483, 'F0025': 0.129}
SAMPLE_ID_CANDIDATES = ['SAMPLE_ID', 'Sample_ID', 'sample_id', 'Sample', 'sample', 'PATIENT_ID', 'Patient_ID', 'patient_id', 'Cell_line']
SUBTYPE_CANDIDATES = ['Molecular_Subtype', 'Molecular Subtype', 'molecular_subtype', 'SUBTYPE', 'Subtype', 'subtype', 'TCGA_Molecular_Subtype', 'TCGA Subtype']
CANONICAL_SUBTYPES = ['POLE-ultramutated', 'MSI-hypermutated', 'Copy-number-low', 'Copy-number-high']

def detect_column(df, candidates, purpose):
    lower_lookup = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in lower_lookup:
            return lower_lookup[key]
    raise ValueError(f'Could not detect {purpose} column.\nAvailable columns include: {list(df.columns[:30])}\nEdit the candidate list in the script if needed.')

def normalize_sample_id(x):
    if pd.isna(x):
        return np.nan
    return str(x).strip().upper()

def normalize_subtype(value):
    if pd.isna(value):
        return np.nan
    original = str(value).strip()
    s = re.sub('[_\\-]+', ' ', original.lower())
    s = re.sub('\\s+', ' ', s).strip()
    if 'pole' in s:
        return 'POLE-ultramutated'
    if 'msi' in s or 'microsatellite' in s:
        return 'MSI-hypermutated'
    if 'copy' in s and 'number' in s and ('high' in s) or 'cn high' in s or 'serous' in s:
        return 'Copy-number-high'
    if 'copy' in s and 'number' in s and ('low' in s) or 'cn low' in s:
        return 'Copy-number-low'
    return original

def read_table(path):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {'.tsv', '.txt'}:
        return pd.read_csv(path, sep='\t', low_memory=False)
    return pd.read_csv(path, low_memory=False)

def prepare_omics_matrix(df, sample_col):
    x = df.copy()
    x[sample_col] = x[sample_col].map(normalize_sample_id)
    x = x.dropna(subset=[sample_col])
    value_cols = [c for c in x.columns if c != sample_col]
    x[value_cols] = x[value_cols].apply(pd.to_numeric, errors='coerce')
    x = x.groupby(sample_col, as_index=True)[value_cols].mean()
    x.columns = [str(c).strip().upper() for c in x.columns]
    if x.columns.duplicated().any():
        x = x.T.groupby(level=0).mean().T
    return x

def zscore_columns(df):
    means = df.mean(axis=0)
    stds = df.std(axis=0, ddof=0).replace(0, np.nan)
    return (df - means) / stds

def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out

def main(args):
    project_root = Path(args.project_root).resolve()
    rna_file = Path(args.rna) if args.rna else project_root / 'data' / 'ucec' / 'tcga_ucec_rna.csv'
    cnv_file = Path(args.cnv) if args.cnv else project_root / 'data' / 'ucec' / 'tcga_ucec_cnv.csv'
    subtype_file = Path(args.subtypes) if args.subtypes else project_root / 'data' / 'ucec' / 'tcga_ucec_subtypes.csv'
    gene_file = Path(args.family_genes) if args.family_genes else project_root / 'results' / 'BRCA_PRISM_step7_family_gene_support.csv'
    results_dir = project_root / 'results'
    figures_dir = project_root / 'figures'
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    print('=' * 72)
    print('BRCA-PRISM: TCGA-UCEC EXPLORATORY CROSS-DISEASE PILOT')
    print('=' * 72)
    print('RNA:', rna_file)
    print('CNV:', cnv_file)
    print('Subtypes:', subtype_file)
    print('Family genes:', gene_file)
    for p in [rna_file, cnv_file, subtype_file, gene_file]:
        if not p.exists():
            raise FileNotFoundError(f'\nRequired input file was not found:\n{p}\n\nSupply the correct path with command-line arguments. Run with --help for details.')
    rna_raw = read_table(rna_file)
    cnv_raw = read_table(cnv_file)
    subtype_raw = read_table(subtype_file)
    family_genes = read_table(gene_file)
    rna_id = detect_column(rna_raw, SAMPLE_ID_CANDIDATES, 'RNA sample-ID')
    cnv_id = detect_column(cnv_raw, SAMPLE_ID_CANDIDATES, 'CNV sample-ID')
    subtype_id = detect_column(subtype_raw, SAMPLE_ID_CANDIDATES, 'subtype sample-ID')
    subtype_col = detect_column(subtype_raw, SUBTYPE_CANDIDATES, 'molecular-subtype')
    labels = subtype_raw[[subtype_id, subtype_col]].copy()
    labels.columns = ['Sample_ID', 'UCEC_Subtype']
    labels['Sample_ID'] = labels['Sample_ID'].map(normalize_sample_id)
    labels['UCEC_Subtype'] = labels['UCEC_Subtype'].map(normalize_subtype)
    labels = labels.dropna().drop_duplicates('Sample_ID')
    print('\nSubtype labels found:')
    print(labels['UCEC_Subtype'].value_counts(dropna=False))
    unexpected = sorted(set(labels['UCEC_Subtype'].dropna().unique()) - set(CANONICAL_SUBTYPES))
    if unexpected:
        print('\nWARNING: Unrecognized subtype labels:')
        for x in unexpected:
            print('  ', x)
        print('These samples will be excluded unless you extend normalize_subtype().')
    labels = labels[labels['UCEC_Subtype'].isin(CANONICAL_SUBTYPES)].copy()
    rna = prepare_omics_matrix(rna_raw, rna_id)
    cnv = prepare_omics_matrix(cnv_raw, cnv_id)
    common_samples = set(rna.index) & set(cnv.index) & set(labels['Sample_ID'])
    common_samples = sorted(common_samples)
    if not common_samples:
        raise ValueError('No common samples were found across RNA, CNV, and subtype files.')
    rna = rna.loc[common_samples]
    cnv = cnv.loc[common_samples]
    labels = labels.set_index('Sample_ID').loc[common_samples].reset_index()
    print('\nAligned UCEC cohort:', len(common_samples), 'samples')
    print(labels['UCEC_Subtype'].value_counts())
    if len(common_samples) != 507:
        print(f'\nNOTE: The paper reports n = 507 tumors with molecular-subtype labels after alignment. Your current files produce n = {len(common_samples)}. Check the cBioPortal export/version if you intend to reproduce the paper exactly.')
    required_cols = {'Family_ID', 'Gene_Symbol', 'Resource_Support'}
    missing = required_cols - set(family_genes.columns)
    if missing:
        raise ValueError('Family-gene file is missing required columns: ' + ', '.join(sorted(missing)))
    family_genes['Family_ID'] = family_genes['Family_ID'].astype(str).str.strip()
    family_genes['Gene_Symbol'] = family_genes['Gene_Symbol'].astype(str).str.upper().str.strip()
    family_genes['Resource_Support'] = pd.to_numeric(family_genes['Resource_Support'], errors='coerce')
    consensus = family_genes[family_genes['Family_ID'].isin(FAMILIES) & (family_genes['Resource_Support'] >= 2)].copy()
    if consensus.empty:
        raise ValueError('No consensus genes were found for F0006/F0025 (Resource_Support >= 2).')
    print('\nConsensus genes by family:')
    print(consensus.groupby('Family_ID')['Gene_Symbol'].nunique())
    rna_z = zscore_columns(rna)
    cnv_z = zscore_columns(cnv)
    activity = labels.copy()
    coverage_rows = []
    for family_id, family_name in FAMILIES.items():
        genes = sorted(consensus.loc[consensus['Family_ID'] == family_id, 'Gene_Symbol'].dropna().unique())
        rna_genes = [g for g in genes if g in rna_z.columns]
        cnv_genes = [g for g in genes if g in cnv_z.columns]
        if not rna_genes and (not cnv_genes):
            raise ValueError(f'No measurable RNA/CNV genes were found for {family_id}.')
        rna_score = rna_z[rna_genes].mean(axis=1) if rna_genes else pd.Series(np.nan, index=rna_z.index)
        cnv_score = cnv_z[cnv_genes].mean(axis=1) if cnv_genes else pd.Series(np.nan, index=cnv_z.index)
        combined = pd.concat([rna_score.rename('RNA'), cnv_score.rename('CNV')], axis=1).mean(axis=1)
        score_map = combined.to_dict()
        activity[f'{family_id}_Combined'] = activity['Sample_ID'].map(score_map)
        coverage_rows.append({'Family_ID': family_id, 'Family_Name': family_name, 'Consensus_Genes': len(genes), 'RNA_Genes_Measured': len(rna_genes), 'CNV_Genes_Measured': len(cnv_genes)})
    coverage = pd.DataFrame(coverage_rows)
    association_rows = []
    for family_id, family_name in FAMILIES.items():
        score_col = f'{family_id}_Combined'
        groups = []
        medians = {}
        for subtype in CANONICAL_SUBTYPES:
            values = activity.loc[activity['UCEC_Subtype'] == subtype, score_col].dropna().to_numpy()
            if len(values) == 0:
                raise ValueError(f'No scores are available for subtype: {subtype}')
            groups.append(values)
            medians[subtype] = float(np.median(values))
        H, p_value = kruskal(*groups)
        n = sum((len(g) for g in groups))
        k = len(groups)
        epsilon2 = max(0.0, (H - k + 1) / (n - k))
        row = {'Family_ID': family_id, 'Family_Name': family_name, 'N': n, 'Kruskal_H': H, 'P_Value': p_value, 'Epsilon_Squared': epsilon2}
        for subtype in CANONICAL_SUBTYPES:
            row[f'Median_{subtype}'] = medians[subtype]
        association_rows.append(row)
    association = pd.DataFrame(association_rows)
    association['FDR_BH_two_family'] = bh_fdr(association['P_Value'])
    comparison = association[['Family_ID', 'Family_Name', 'Epsilon_Squared', 'P_Value']].copy()
    comparison = comparison.rename(columns={'Epsilon_Squared': 'TCGA_UCEC_Epsilon2', 'P_Value': 'TCGA_UCEC_P_Value'})
    comparison['TCGA_BRCA_Epsilon2'] = comparison['Family_ID'].map(BRCA_EFFECT_SIZE)
    comparison['Absolute_Change'] = comparison['TCGA_UCEC_Epsilon2'] - comparison['TCGA_BRCA_Epsilon2']
    comparison['Percent_Change'] = 100 * comparison['Absolute_Change'] / comparison['TCGA_BRCA_Epsilon2']
    median_rows = []
    for family_id, family_name in FAMILIES.items():
        score_col = f'{family_id}_Combined'
        for subtype in CANONICAL_SUBTYPES:
            vals = activity.loc[activity['UCEC_Subtype'] == subtype, score_col]
            median_rows.append({'Family_ID': family_id, 'Family_Name': family_name, 'UCEC_Subtype': subtype, 'Median_Activity': vals.median(), 'N': vals.notna().sum()})
    medians = pd.DataFrame(median_rows)
    activity_file = results_dir / 'BRCA_PRISM_UCEC_two_family_activity_scores.csv'
    association_file = results_dir / 'BRCA_PRISM_UCEC_two_family_association.csv'
    comparison_file = results_dir / 'BRCA_PRISM_UCEC_vs_BRCA_effect_size_comparison.csv'
    median_file = results_dir / 'BRCA_PRISM_UCEC_subtype_median_activity.csv'
    coverage_file = results_dir / 'BRCA_PRISM_UCEC_gene_coverage.csv'
    activity.to_csv(activity_file, index=False)
    association.to_csv(association_file, index=False)
    comparison.to_csv(comparison_file, index=False)
    medians.to_csv(median_file, index=False)
    coverage.to_csv(coverage_file, index=False)
    x = np.arange(len(comparison))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width / 2, comparison['TCGA_BRCA_Epsilon2'], width, label='TCGA-BRCA')
    ax.bar(x + width / 2, comparison['TCGA_UCEC_Epsilon2'], width, label='TCGA-UCEC')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{fid}\n{name}' for fid, name in zip(comparison['Family_ID'], comparison['Family_Name'])])
    ax.set_ylabel('Kruskal-Wallis epsilon-squared')
    ax.set_title('Exploratory cross-disease pathway-family comparison')
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / 'BRCA_PRISM_UCEC_two_family_effect_sizes.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    pivot = medians.pivot(index='Family_ID', columns='UCEC_Subtype', values='Median_Activity')
    pivot = pivot.reindex(index=list(FAMILIES.keys()), columns=CANONICAL_SUBTYPES)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(pivot.to_numpy(), aspect='auto')
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([f'{fid} — {FAMILIES[fid]}' for fid in pivot.index])
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=25, ha='right')
    ax.set_title('TCGA-UCEC median consensus-family activity by subtype')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Median standardized family activity')
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            if pd.notna(value):
                ax.text(j, i, f'{value:.2f}', ha='center', va='center')
    fig.tight_layout()
    fig.savefig(figures_dir / 'BRCA_PRISM_UCEC_subtype_activity.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('\n' + '=' * 72)
    print('TCGA-UCEC RESULTS')
    print('=' * 72)
    print(association[['Family_ID', 'Family_Name', 'N', 'Kruskal_H', 'P_Value', 'Epsilon_Squared']].to_string(index=False))
    print('\nTCGA-BRCA vs TCGA-UCEC:')
    print(comparison.to_string(index=False))
    print('F0006: UCEC epsilon² ≈ 0.164, H ≈ 85.70, p ≈ 1.8e-18')
    print('F0025: UCEC epsilon² ≈ 0.115, H ≈ 60.98, p ≈ 3.6e-13')
    print('\nSaved:')
    for p in [activity_file, association_file, comparison_file, median_file, coverage_file, figures_dir / 'BRCA_PRISM_UCEC_two_family_effect_sizes.png', figures_dir / 'BRCA_PRISM_UCEC_subtype_activity.png']:
        print(' ', p)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BRCA-PRISM exploratory TCGA-UCEC cross-disease analysis for F0006 and F0025.')
    parser.add_argument('--project-root', default='.', help='Repository/project root. Default: current directory.')
    parser.add_argument('--rna', default=None, help='Path to TCGA-UCEC RNA sample-by-gene table.')
    parser.add_argument('--cnv', default=None, help='Path to TCGA-UCEC CNV sample-by-gene table.')
    parser.add_argument('--subtypes', default=None, help='Path to TCGA-UCEC molecular-subtype table.')
    parser.add_argument('--family-genes', default=None, help='Path to BRCA_PRISM_step7_family_gene_support.csv. Default: results/BRCA_PRISM_step7_family_gene_support.csv')
    main(parser.parse_args())
