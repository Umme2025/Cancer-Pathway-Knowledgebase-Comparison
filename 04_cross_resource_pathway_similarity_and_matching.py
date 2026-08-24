pathway_gene_sets = {}
for (db, pid, pname), group in master_df.groupby(['Database', 'Pathway_ID', 'Pathway_Name']):
    pathway_gene_sets[db, pid, pname] = set(group['Gene_Symbol'].dropna())
print('Total pathway gene sets:', len(pathway_gene_sets))

def pathway_similarity(reference_genes, candidate_genes):
    shared = reference_genes & candidate_genes
    n_shared = len(shared)
    n_ref = len(reference_genes)
    n_candidate = len(candidate_genes)
    precision = n_shared / n_candidate if n_candidate > 0 else 0
    recall = n_shared / n_ref if n_ref > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0
    union = reference_genes | candidate_genes
    jaccard = n_shared / len(union) if len(union) > 0 else 0
    overlap_coefficient = n_shared / min(n_ref, n_candidate) if min(n_ref, n_candidate) > 0 else 0
    return {'Shared_Genes': n_shared, 'Precision': precision, 'Recall': recall, 'F1': f1, 'Jaccard': jaccard, 'Overlap_Coefficient': overlap_coefficient}
reference_db = 'LCPathways'
target_dbs = ['KEGG', 'Reactome', 'WikiPathways', 'MSigDB_Hallmark']
lc_pathways = [key for key in pathway_gene_sets if key[0] == reference_db]
best_match_rows = []
for ref_db, ref_id, ref_name in lc_pathways:
    ref_genes = pathway_gene_sets[ref_db, ref_id, ref_name]
    for target_db in target_dbs:
        best_result = None
        for (db, pid, pname), candidate_genes in pathway_gene_sets.items():
            if db != target_db:
                continue
            metrics = pathway_similarity(ref_genes, candidate_genes)
            result = {'Reference_Database': ref_db, 'Reference_Pathway_ID': ref_id, 'Reference_Pathway_Name': ref_name, 'Reference_Gene_Count': len(ref_genes), 'Target_Database': db, 'Best_Match_Pathway_ID': pid, 'Best_Match_Pathway_Name': pname, 'Target_Gene_Count': len(candidate_genes), **metrics}
            if best_result is None or result['F1'] > best_result['F1']:
                best_result = result
        best_match_rows.append(best_result)
best_match_df = pd.DataFrame(best_match_rows)
print(best_match_df.shape)
print(best_match_df.head())
metric_cols = ['Precision', 'Recall', 'F1', 'Jaccard', 'Overlap_Coefficient']
best_match_df[metric_cols] = best_match_df[metric_cols].round(4)
best_match_df.to_csv(Path.cwd() / 'results' / 'LCPathways_best_matching_pathways_across_resources.csv', index=False)
print('Saved.')
best_match_summary = best_match_df.groupby('Target_Database').agg(Mean_F1=('F1', 'mean'), Median_F1=('F1', 'median'), Mean_Jaccard=('Jaccard', 'mean'), Median_Jaccard=('Jaccard', 'median'), Mean_Recall=('Recall', 'mean'), Mean_Precision=('Precision', 'mean'), Mean_Overlap_Coefficient=('Overlap_Coefficient', 'mean')).round(4)
print(best_match_summary)
import matplotlib.pyplot as plt
order = ['KEGG', 'Reactome', 'WikiPathways', 'MSigDB_Hallmark']
data = [best_match_df.loc[best_match_df['Target_Database'] == db, 'F1'].values for db in order]
plt.figure(figsize=(9, 6))
plt.boxplot(data, tick_labels=order, showfliers=False)
plt.ylabel('Best-match F1 score')
plt.xlabel('Target pathway resource')
plt.title('LCPathways best-match similarity across pathway resources')
plt.tight_layout()
plt.savefig(Path.cwd() / 'figures' / 'LCPathways_best_match_F1_distribution.png', dpi=300, bbox_inches='tight')
plt.show()

def classify_f1(x):
    if x >= 0.5:
        return 'Strong (F1 >= 0.50)'
    elif x >= 0.25:
        return 'Moderate (0.25 <= F1 < 0.50)'
    else:
        return 'Weak (F1 < 0.25)'
best_match_df['Match_Strength'] = best_match_df['F1'].apply(classify_f1)
match_strength_summary = best_match_df.groupby(['Target_Database', 'Match_Strength']).size().reset_index(name='Pathway_Count')
print(match_strength_summary)
top_matches = best_match_df.sort_values(['Target_Database', 'F1'], ascending=[True, False]).groupby('Target_Database').head(10)
display(top_matches[['Reference_Pathway_Name', 'Target_Database', 'Best_Match_Pathway_Name', 'Reference_Gene_Count', 'Target_Gene_Count', 'Shared_Genes', 'Precision', 'Recall', 'F1', 'Jaccard']])
top_matches.to_csv(Path.cwd() / 'results' / 'LCPathways_top10_best_matches_per_resource.csv', index=False)
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
RESULTS_DIR = PROJECT_ROOT / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
master = pd.read_csv(MASTER_FILE)
resource_order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB Hallmark']
master = master.dropna(subset=['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']).copy()
master['Database'] = master['Database'].astype(str).str.strip()
master['Pathway_ID'] = master['Pathway_ID'].astype(str).str.strip()
master['Pathway_Name'] = master['Pathway_Name'].astype(str).str.strip()
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master = master[master['Database'].isin(resource_order)].copy()
master = master.drop_duplicates(subset=['Database', 'Pathway_ID', 'Gene_Symbol'])
print('Master shape:', master.shape)
print('Resources:', master['Database'].unique().tolist())
resource_pathways = {}
for resource in resource_order:
    df = master[master['Database'] == resource].copy()
    pathway_info = df[['Pathway_ID', 'Pathway_Name']].drop_duplicates(subset=['Pathway_ID']).reset_index(drop=True)
    pathway_gene_sets = df.groupby('Pathway_ID')['Gene_Symbol'].apply(set).to_dict()
    resource_pathways[resource] = {'info': pathway_info, 'gene_sets': pathway_gene_sets}
    print(f'{resource:<18}: {len(pathway_info):,} pathways')

def directional_best_matches(source_resource, target_resource, source_data, target_data):
    print('\n' + '=' * 90)
    print(f'{source_resource}  -->  {target_resource}')
    print('=' * 90)
    source_info = source_data['info']
    target_info = target_data['info']
    source_gene_sets = source_data['gene_sets']
    target_gene_sets = target_data['gene_sets']
    all_genes = set()
    for genes in source_gene_sets.values():
        all_genes.update(genes)
    for genes in target_gene_sets.values():
        all_genes.update(genes)
    all_genes = sorted(all_genes)
    gene_index = {gene: i for i, gene in enumerate(all_genes)}
    source_ids = source_info['Pathway_ID'].tolist()
    source_index = {pathway_id: i for i, pathway_id in enumerate(source_ids)}
    src_rows = []
    src_cols = []
    for pathway_id, genes in source_gene_sets.items():
        i = source_index[pathway_id]
        for gene in genes:
            src_rows.append(i)
            src_cols.append(gene_index[gene])
    X_source = csr_matrix((np.ones(len(src_rows), dtype=np.int32), (src_rows, src_cols)), shape=(len(source_ids), len(all_genes)))
    X_source.data[:] = 1
    target_ids = target_info['Pathway_ID'].tolist()
    target_index = {pathway_id: i for i, pathway_id in enumerate(target_ids)}
    tgt_rows = []
    tgt_cols = []
    for pathway_id, genes in target_gene_sets.items():
        j = target_index[pathway_id]
        for gene in genes:
            tgt_rows.append(j)
            tgt_cols.append(gene_index[gene])
    X_target = csr_matrix((np.ones(len(tgt_rows), dtype=np.int32), (tgt_rows, tgt_cols)), shape=(len(target_ids), len(all_genes)))
    X_target.data[:] = 1
    source_sizes = np.asarray(X_source.sum(axis=1)).ravel()
    target_sizes = np.asarray(X_target.sum(axis=1)).ravel()
    shared_matrix = (X_source @ X_target.T).tocsr()
    results = []
    for i, source_id in enumerate(source_ids):
        start = shared_matrix.indptr[i]
        end = shared_matrix.indptr[i + 1]
        candidate_indices = shared_matrix.indices[start:end]
        shared_counts = shared_matrix.data[start:end]
        best_idx = None
        best_shared = 0
        best_jaccard = 0.0
        best_overlap = 0.0
        if len(candidate_indices) == 0:
            source_name = source_info.loc[i, 'Pathway_Name']
            results.append({'Source_Resource': source_resource, 'Target_Resource': target_resource, 'Source_Pathway_ID': source_id, 'Source_Pathway_Name': source_name, 'Source_Pathway_Size': int(source_sizes[i]), 'Target_Pathway_ID': None, 'Target_Pathway_Name': None, 'Target_Pathway_Size': np.nan, 'Shared_Genes': 0, 'Jaccard': 0.0, 'Overlap_Coefficient': 0.0})
            continue
        for j, shared in zip(candidate_indices, shared_counts):
            shared = int(shared)
            source_size = int(source_sizes[i])
            target_size = int(target_sizes[j])
            union_size = source_size + target_size - shared
            min_size = min(source_size, target_size)
            jaccard = shared / union_size if union_size > 0 else 0.0
            overlap = shared / min_size if min_size > 0 else 0.0
            better = False
            if jaccard > best_jaccard:
                better = True
            elif np.isclose(jaccard, best_jaccard):
                if overlap > best_overlap:
                    better = True
                elif np.isclose(overlap, best_overlap):
                    if shared > best_shared:
                        better = True
            if better:
                best_idx = j
                best_shared = shared
                best_jaccard = jaccard
                best_overlap = overlap
        source_name = source_info.loc[i, 'Pathway_Name']
        if best_idx is None:
            results.append({'Source_Resource': source_resource, 'Target_Resource': target_resource, 'Source_Pathway_ID': source_id, 'Source_Pathway_Name': source_name, 'Source_Pathway_Size': int(source_sizes[i]), 'Target_Pathway_ID': None, 'Target_Pathway_Name': None, 'Target_Pathway_Size': np.nan, 'Shared_Genes': 0, 'Jaccard': 0.0, 'Overlap_Coefficient': 0.0})
        else:
            target_name = target_info.loc[best_idx, 'Pathway_Name']
            results.append({'Source_Resource': source_resource, 'Target_Resource': target_resource, 'Source_Pathway_ID': source_id, 'Source_Pathway_Name': source_name, 'Source_Pathway_Size': int(source_sizes[i]), 'Target_Pathway_ID': target_ids[best_idx], 'Target_Pathway_Name': target_name, 'Target_Pathway_Size': int(target_sizes[best_idx]), 'Shared_Genes': int(best_shared), 'Jaccard': float(best_jaccard), 'Overlap_Coefficient': float(best_overlap)})
    result_df = pd.DataFrame(results)
    print('Source pathways:', len(result_df))
    print('Mean best Jaccard:', round(result_df['Jaccard'].mean(), 4))
    print('Median best Jaccard:', round(result_df['Jaccard'].median(), 4))
    print('Mean best overlap coefficient:', round(result_df['Overlap_Coefficient'].mean(), 4))
    return result_df
all_directional_results = []
for source_resource in resource_order:
    for target_resource in resource_order:
        if source_resource == target_resource:
            continue
        result_df = directional_best_matches(source_resource, target_resource, resource_pathways[source_resource], resource_pathways[target_resource])
        all_directional_results.append(result_df)
mapping_df = pd.concat(all_directional_results, ignore_index=True)
print('\n' + '=' * 100)
print('ALL DIRECTIONAL BEST-MATCH RESULTS')
print('=' * 100)
print('Total source->target pathway mappings:', len(mapping_df))
display(mapping_df.head(20))
directional_summary = mapping_df.groupby(['Source_Resource', 'Target_Resource']).agg(Number_of_Source_Pathways=('Source_Pathway_ID', 'count'), Mean_Best_Jaccard=('Jaccard', 'mean'), Median_Best_Jaccard=('Jaccard', 'median'), Mean_Best_Overlap=('Overlap_Coefficient', 'mean'), Median_Best_Overlap=('Overlap_Coefficient', 'median'), Fraction_Jaccard_ge_0_25=('Jaccard', lambda x: (x >= 0.25).mean()), Fraction_Jaccard_ge_0_50=('Jaccard', lambda x: (x >= 0.5).mean()), Fraction_Overlap_ge_0_50=('Overlap_Coefficient', lambda x: (x >= 0.5).mean()), Fraction_Overlap_ge_0_75=('Overlap_Coefficient', lambda x: (x >= 0.75).mean())).reset_index()
print('\n' + '=' * 100)
print('DIRECTIONAL PATHWAY-SIMILARITY SUMMARY')
print('=' * 100)
display(directional_summary.round(4))
mean_jaccard_matrix = pd.DataFrame(np.nan, index=resource_order, columns=resource_order)
median_jaccard_matrix = pd.DataFrame(np.nan, index=resource_order, columns=resource_order)
mean_overlap_matrix = pd.DataFrame(np.nan, index=resource_order, columns=resource_order)
for _, row in directional_summary.iterrows():
    source = row['Source_Resource']
    target = row['Target_Resource']
    mean_jaccard_matrix.loc[source, target] = row['Mean_Best_Jaccard']
    median_jaccard_matrix.loc[source, target] = row['Median_Best_Jaccard']
    mean_overlap_matrix.loc[source, target] = row['Mean_Best_Overlap']
print('\nMean best-Jaccard matrix:')
display(mean_jaccard_matrix.round(4))
print('\nMedian best-Jaccard matrix:')
display(median_jaccard_matrix.round(4))
print('\nMean best-overlap-coefficient matrix:')
display(mean_overlap_matrix.round(4))
top_candidates = mapping_df.sort_values(['Jaccard', 'Overlap_Coefficient', 'Shared_Genes'], ascending=[False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('TOP 50 CROSS-RESOURCE CANDIDATE MATCHES')
print('=' * 100)
display(top_candidates[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient']].head(50).round(4))
mapping_file = RESULTS_DIR / 'BRCA_PRISM_step4A_directional_best_pathway_matches.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_step4A_directional_similarity_summary.csv'
mean_jaccard_file = RESULTS_DIR / 'BRCA_PRISM_step4A_mean_best_jaccard_matrix.csv'
median_jaccard_file = RESULTS_DIR / 'BRCA_PRISM_step4A_median_best_jaccard_matrix.csv'
mean_overlap_file = RESULTS_DIR / 'BRCA_PRISM_step4A_mean_best_overlap_matrix.csv'
top_candidates_file = RESULTS_DIR / 'BRCA_PRISM_step4A_top_candidate_matches.csv'
mapping_df.to_csv(mapping_file, index=False)
directional_summary.to_csv(summary_file, index=False)
mean_jaccard_matrix.to_csv(mean_jaccard_file)
median_jaccard_matrix.to_csv(median_jaccard_file)
mean_overlap_matrix.to_csv(mean_overlap_file)
top_candidates.to_csv(top_candidates_file, index=False)
print('\n' + '=' * 100)
print('STEP 4A COMPLETE')
print('=' * 100)
print('\nSaved:')
print(mapping_file)
print(summary_file)
print(mean_jaccard_file)
print(median_jaccard_file)
print(mean_overlap_file)
print(top_candidates_file)
print('\nPlease send me:')
print('1. DIRECTIONAL PATHWAY-SIMILARITY SUMMARY')
print('2. Mean best-Jaccard matrix')
print('3. Mean best-overlap-coefficient matrix')
print('4. First ~20 rows of TOP candidate matches')
from pathlib import Path
import pandas as pd
import numpy as np
import re
import unicodedata
from rapidfuzz.fuzz import ratio, token_sort_ratio
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
INPUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step4A_directional_best_pathway_matches.csv'
print('Input file:')
print(INPUT_FILE)
mapping = pd.read_csv(INPUT_FILE)
print('\nStep 4A mapping shape:', mapping.shape)
required_columns = ['Source_Resource', 'Target_Resource', 'Source_Pathway_ID', 'Source_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_ID', 'Target_Pathway_Name', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient']
missing = [c for c in required_columns if c not in mapping.columns]
if missing:
    raise ValueError(f'Missing required columns: {missing}')

def normalize_pathway_name(name):
    if pd.isna(name):
        return ''
    name = str(name)
    name = unicodedata.normalize('NFKD', name)
    name = name.lower()
    name = name.replace('&', ' and ')
    name = re.sub('[_\\-/,:;()\\[\\]{}]+', ' ', name)
    name = re.sub('[^a-z0-9\\s]', ' ', name)
    name = re.sub('\\s+', ' ', name).strip()
    return name
mapping['Source_Name_Normalized'] = mapping['Source_Pathway_Name'].apply(normalize_pathway_name)
mapping['Target_Name_Normalized'] = mapping['Target_Pathway_Name'].apply(normalize_pathway_name)

def lexical_similarity(row):
    a = row['Source_Name_Normalized']
    b = row['Target_Name_Normalized']
    if not a or not b:
        return 0.0
    return ratio(a, b) / 100.0
mapping['Name_Levenshtein_Similarity'] = mapping.apply(lexical_similarity, axis=1)

def token_similarity(row):
    a = row['Source_Name_Normalized']
    b = row['Target_Name_Normalized']
    if not a or not b:
        return 0.0
    return token_sort_ratio(a, b) / 100.0
mapping['Name_Token_Similarity'] = mapping.apply(token_similarity, axis=1)
mapping['Exact_Normalized_Name_Match'] = mapping['Source_Name_Normalized'] == mapping['Target_Name_Normalized']
mapping['Name_Similarity'] = mapping[['Name_Levenshtein_Similarity', 'Name_Token_Similarity']].max(axis=1)
forward_keys = set(zip(mapping['Source_Resource'], mapping['Source_Pathway_ID'].astype(str), mapping['Target_Resource'], mapping['Target_Pathway_ID'].astype(str)))
reciprocal_flags = []
for _, row in mapping.iterrows():
    reverse_key = (row['Target_Resource'], str(row['Target_Pathway_ID']), row['Source_Resource'], str(row['Source_Pathway_ID']))
    reciprocal_flags.append(reverse_key in forward_keys)
mapping['Reciprocal_Best_Match'] = reciprocal_flags
mapping['High_Gene_Similarity'] = mapping['Jaccard'] >= 0.5
mapping['High_Containment'] = mapping['Overlap_Coefficient'] >= 0.75
mapping['High_Name_Similarity'] = mapping['Name_Similarity'] >= 0.8

def evidence_level(row):
    gene_good = row['Jaccard'] >= 0.5 or row['Overlap_Coefficient'] >= 0.75
    name_good = row['Name_Similarity'] >= 0.8
    reciprocal = row['Reciprocal_Best_Match']
    if reciprocal and gene_good and name_good:
        return 'Strong candidate'
    elif gene_good and name_good or (reciprocal and gene_good) or (reciprocal and name_good):
        return 'Moderate candidate'
    else:
        return 'Weak candidate'
mapping['Candidate_Evidence'] = mapping.apply(evidence_level, axis=1)
print('\n' + '=' * 100)
print('STEP 4B NAME-SIMILARITY SUMMARY')
print('=' * 100)
print('Mean name similarity:', round(mapping['Name_Similarity'].mean(), 4))
print('Median name similarity:', round(mapping['Name_Similarity'].median(), 4))
print('Exact normalized-name matches:', int(mapping['Exact_Normalized_Name_Match'].sum()))
print('Reciprocal best matches:', int(mapping['Reciprocal_Best_Match'].sum()))
print('Fraction reciprocal:', round(mapping['Reciprocal_Best_Match'].mean(), 4))
evidence_summary = mapping['Candidate_Evidence'].value_counts().rename_axis('Candidate_Evidence').reset_index(name='Number_of_Mappings')
evidence_summary['Percentage'] = evidence_summary['Number_of_Mappings'] / len(mapping) * 100
print('\nCandidate evidence summary:')
display(evidence_summary.round(2))
pair_summary = mapping.groupby(['Source_Resource', 'Target_Resource']).agg(Number_of_Source_Pathways=('Source_Pathway_ID', 'count'), Mean_Jaccard=('Jaccard', 'mean'), Mean_Overlap=('Overlap_Coefficient', 'mean'), Mean_Name_Similarity=('Name_Similarity', 'mean'), Fraction_NameSimilarity_ge_0_80=('Name_Similarity', lambda x: (x >= 0.8).mean()), Fraction_Reciprocal=('Reciprocal_Best_Match', 'mean'), Fraction_Strong_Candidates=('Candidate_Evidence', lambda x: (x == 'Strong candidate').mean())).reset_index()
print('\n' + '=' * 100)
print('RESOURCE-PAIR NAME + GENE SIMILARITY SUMMARY')
print('=' * 100)
display(pair_summary.round(4))
strong_candidates = mapping[mapping['Candidate_Evidence'] == 'Strong candidate'].sort_values(['Jaccard', 'Name_Similarity', 'Overlap_Coefficient', 'Shared_Genes'], ascending=[False, False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('TOP STRONG CANDIDATE MAPPINGS')
print('=' * 100)
display(strong_candidates[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']].head(50).round(4))
suspicious = mapping[(mapping['Jaccard'] >= 0.4) & (mapping['Name_Similarity'] < 0.5)].sort_values('Jaccard', ascending=False).reset_index(drop=True)
print('\n' + '=' * 100)
print('HIGH GENE OVERLAP BUT LOW NAME SIMILARITY')
print('=' * 100)
display(suspicious[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']].head(30).round(4))
name_gene_disagreement = mapping[(mapping['Name_Similarity'] >= 0.8) & (mapping['Jaccard'] < 0.25)].sort_values(['Name_Similarity', 'Jaccard'], ascending=[False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('HIGH NAME SIMILARITY BUT LOW GENE JACCARD')
print('=' * 100)
display(name_gene_disagreement[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']].head(30).round(4))
full_file = RESULTS_DIR / 'BRCA_PRISM_step4B_pathway_matches_with_name_similarity.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_step4B_resource_pair_summary.csv'
strong_file = RESULTS_DIR / 'BRCA_PRISM_step4B_strong_candidate_mappings.csv'
suspicious_file = RESULTS_DIR / 'BRCA_PRISM_step4B_high_gene_low_name_candidates.csv'
name_disagreement_file = RESULTS_DIR / 'BRCA_PRISM_step4B_high_name_low_gene_candidates.csv'
evidence_file = RESULTS_DIR / 'BRCA_PRISM_step4B_candidate_evidence_summary.csv'
mapping.to_csv(full_file, index=False)
pair_summary.to_csv(summary_file, index=False)
strong_candidates.to_csv(strong_file, index=False)
suspicious.to_csv(suspicious_file, index=False)
name_gene_disagreement.to_csv(name_disagreement_file, index=False)
evidence_summary.to_csv(evidence_file, index=False)
print('\n' + '=' * 100)
print('STEP 4B COMPLETE')
print('=' * 100)
print('\nSaved:')
print(full_file)
print(summary_file)
print(strong_file)
print(suspicious_file)
print(name_disagreement_file)
print(evidence_file)
from pathlib import Path
import pandas as pd
import numpy as np
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
INPUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step4B_pathway_matches_with_name_similarity.csv'
print('Input:')
print(INPUT_FILE)
df = pd.read_csv(INPUT_FILE)
print('\nShape:', df.shape)
required_columns = ['Source_Resource', 'Target_Resource', 'Source_Pathway_ID', 'Source_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_ID', 'Target_Pathway_Name', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']
missing = [c for c in required_columns if c not in df.columns]
if missing:
    raise ValueError(f'Missing required columns: {missing}')
df['Smaller_Pathway_Size'] = df[['Source_Pathway_Size', 'Target_Pathway_Size']].min(axis=1)
df['Larger_Pathway_Size'] = df[['Source_Pathway_Size', 'Target_Pathway_Size']].max(axis=1)
df['Size_Ratio'] = df['Larger_Pathway_Size'] / df['Smaller_Pathway_Size']
df['Absolute_Size_Difference'] = (df['Source_Pathway_Size'] - df['Target_Pathway_Size']).abs()

def containment_direction(row):
    if row['Overlap_Coefficient'] < 0.75:
        return 'No strong containment'
    if row['Source_Pathway_Size'] < row['Target_Pathway_Size']:
        return 'Source may be part of Target'
    elif row['Source_Pathway_Size'] > row['Target_Pathway_Size']:
        return 'Target may be part of Source'
    else:
        return 'Similar size'
df['Containment_Direction'] = df.apply(containment_direction, axis=1)

def classify_relationship(row):
    j = float(row['Jaccard'])
    o = float(row['Overlap_Coefficient'])
    n = float(row['Name_Similarity'])
    reciprocal = bool(row['Reciprocal_Best_Match'])
    shared = int(row['Shared_Genes'])
    source_size = float(row['Source_Pathway_Size'])
    target_size = float(row['Target_Pathway_Size'])
    smaller_size = min(source_size, target_size)
    size_ratio = float(row['Size_Ratio'])
    if reciprocal and shared >= 5 and (j >= 0.5) and (n >= 0.8) and (size_ratio <= 1.5):
        return 'Equivalent candidate'
    if shared >= 3 and smaller_size >= 3 and (o >= 0.75) and (size_ratio >= 1.5) and (n >= 0.4 or reciprocal):
        return 'Hierarchical / isPartOf candidate'
    if shared >= 5 and smaller_size >= 5 and (o >= 0.8) and (size_ratio >= 1.3) and (n >= 0.6) and (j < 0.6):
        return 'Hierarchical / isPartOf candidate'
    if shared >= 3 and (j >= 0.25 and n >= 0.4 or (o >= 0.5 and n >= 0.5) or (reciprocal and j >= 0.2)):
        return 'Related candidate'
    return 'No reliable counterpart'
df['Provisional_Relationship'] = df.apply(classify_relationship, axis=1)
df['High_Gene_Low_Name_Flag'] = (df['Jaccard'] >= 0.4) & (df['Name_Similarity'] < 0.5)
df['High_Name_Low_Gene_Flag'] = (df['Name_Similarity'] >= 0.8) & (df['Jaccard'] < 0.25)
df['Very_High_Containment_Flag'] = df['Overlap_Coefficient'] >= 0.9
df['Tiny_Pathway_Flag'] = df['Smaller_Pathway_Size'] < 3

def review_priority(row):
    relationship = row['Provisional_Relationship']
    if relationship in ['Equivalent candidate', 'Hierarchical / isPartOf candidate']:
        return 'High'
    if row['High_Gene_Low_Name_Flag'] or row['High_Name_Low_Gene_Flag']:
        return 'High'
    if relationship == 'Related candidate':
        return 'Medium'
    return 'Low'
df['Review_Priority'] = df.apply(review_priority, axis=1)
df['Final_Relationship'] = ''
df['Biological_Scope_Compatible'] = ''
df['Biological_Context_Compatible'] = ''
df['Hierarchy_Checked'] = ''
df['Curator_Notes'] = ''
relationship_summary = df['Provisional_Relationship'].value_counts().rename_axis('Provisional_Relationship').reset_index(name='Number_of_Mappings')
relationship_summary['Percentage'] = relationship_summary['Number_of_Mappings'] / len(df) * 100
print('\n' + '=' * 100)
print('REFINED PROVISIONAL RELATIONSHIP SUMMARY')
print('=' * 100)
display(relationship_summary.round(2))
pair_relationship_summary = df.groupby(['Source_Resource', 'Target_Resource', 'Provisional_Relationship']).size().reset_index(name='Number_of_Mappings')
print('\n' + '=' * 100)
print('RELATIONSHIP SUMMARY BY RESOURCE PAIR')
print('=' * 100)
display(pair_relationship_summary)
equivalent_candidates = df[df['Provisional_Relationship'] == 'Equivalent candidate'].sort_values(['Jaccard', 'Name_Similarity', 'Overlap_Coefficient', 'Shared_Genes'], ascending=[False, False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('TOP EQUIVALENT CANDIDATES')
print('=' * 100)
display(equivalent_candidates[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']].head(30).round(4))
hierarchical_candidates = df[df['Provisional_Relationship'] == 'Hierarchical / isPartOf candidate'].sort_values(['Overlap_Coefficient', 'Shared_Genes', 'Name_Similarity', 'Size_Ratio'], ascending=[False, False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('TOP REFINED HIERARCHICAL / isPartOf CANDIDATES')
print('=' * 100)
display(hierarchical_candidates[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Source_Pathway_Size', 'Target_Pathway_Size', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Size_Ratio', 'Containment_Direction']].head(30).round(4))
related_candidates = df[df['Provisional_Relationship'] == 'Related candidate'].sort_values(['Jaccard', 'Overlap_Coefficient', 'Name_Similarity'], ascending=[False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('TOP RELATED CANDIDATES')
print('=' * 100)
display(related_candidates[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match']].head(30).round(4))
manual_review = df[df['Review_Priority'] == 'High'].copy()
manual_review = manual_review.sort_values(['Provisional_Relationship', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity'], ascending=[True, False, False, False]).reset_index(drop=True)
print('\n' + '=' * 100)
print('HIGH-PRIORITY MANUAL REVIEW QUEUE')
print('=' * 100)
print('Mappings requiring high-priority review:', len(manual_review))
display(manual_review[['Source_Resource', 'Target_Resource', 'Source_Pathway_Name', 'Target_Pathway_Name', 'Provisional_Relationship', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Size_Ratio', 'Containment_Direction', 'Tiny_Pathway_Flag', 'High_Gene_Low_Name_Flag', 'High_Name_Low_Gene_Flag']].head(50).round(4))
tiny_summary = df[df['Tiny_Pathway_Flag']]['Provisional_Relationship'].value_counts().rename_axis('Provisional_Relationship').reset_index(name='Number_of_Tiny_Pathway_Mappings')
print('\n' + '=' * 100)
print('TINY-PATHWAY DIAGNOSTIC')
print('=' * 100)
display(tiny_summary)
full_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_provisional_pathway_relationships.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_relationship_summary.csv'
pair_summary_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_relationship_summary_by_resource_pair.csv'
equivalent_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_equivalent_candidates.csv'
hierarchical_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_hierarchical_candidates.csv'
related_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_related_candidates.csv'
manual_review_file = RESULTS_DIR / 'BRCA_PRISM_step5_refined_manual_review_queue.csv'
df.to_csv(full_file, index=False)
relationship_summary.to_csv(summary_file, index=False)
pair_relationship_summary.to_csv(pair_summary_file, index=False)
equivalent_candidates.to_csv(equivalent_file, index=False)
hierarchical_candidates.to_csv(hierarchical_file, index=False)
related_candidates.to_csv(related_file, index=False)
manual_review.to_csv(manual_review_file, index=False)
print('\n' + '=' * 100)
print('STEP 5 REFINED SCREENING COMPLETE')
print('=' * 100)
print('\nSaved:')
print(full_file)
print(summary_file)
print(pair_summary_file)
print(equivalent_file)
print(hierarchical_file)
print(related_file)
print(manual_review_file)
