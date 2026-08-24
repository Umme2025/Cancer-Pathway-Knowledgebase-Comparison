from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / 'data_processed'
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
MASTER_FILE = PROCESSED_DIR / 'BRCA_PRISM_all_pathway_resources_master.csv'
master = pd.read_csv(MASTER_FILE)
required_columns = ['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']
missing = [col for col in required_columns if col not in master.columns]
if missing:
    raise ValueError(f'Missing required columns: {missing}')
print('Master shape:', master.shape)
print('Resources:', master['Database'].unique().tolist())
resource_order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB Hallmark']
pathway_sizes = master.groupby(['Database', 'Pathway_ID', 'Pathway_Name'], as_index=False).agg(Pathway_Size=('Gene_Symbol', 'nunique'))
print('\nIndividual pathway-size table:')
print('Shape:', pathway_sizes.shape)
display(pathway_sizes.head())
structural_summary = pathway_sizes.groupby('Database').agg(Number_of_Pathways=('Pathway_ID', 'nunique'), Mean_Pathway_Size=('Pathway_Size', 'mean'), Median_Pathway_Size=('Pathway_Size', 'median'), Min_Pathway_Size=('Pathway_Size', 'min'), Max_Pathway_Size=('Pathway_Size', 'max'), SD_Pathway_Size=('Pathway_Size', 'std')).reset_index()
unique_gene_counts = master.groupby('Database')['Gene_Symbol'].nunique().rename('Number_of_Unique_Genes').reset_index()
structural_summary = structural_summary.merge(unique_gene_counts, on='Database', how='left')
pathway_size_categories = pathway_sizes.assign(Small_10=lambda x: x['Pathway_Size'] <= 10, Small_25=lambda x: x['Pathway_Size'] <= 25, Medium_26_100=lambda x: (x['Pathway_Size'] >= 26) & (x['Pathway_Size'] <= 100), Large_100=lambda x: x['Pathway_Size'] > 100, Very_Large_500=lambda x: x['Pathway_Size'] > 500).groupby('Database').agg(Fraction_Pathways_10_or_Fewer=('Small_10', 'mean'), Fraction_Pathways_25_or_Fewer=('Small_25', 'mean'), Fraction_Pathways_26_to_100=('Medium_26_100', 'mean'), Fraction_Pathways_Over_100=('Large_100', 'mean'), Fraction_Pathways_Over_500=('Very_Large_500', 'mean')).reset_index()
structural_summary = structural_summary.merge(pathway_size_categories, on='Database', how='left')
gene_depth = master.groupby(['Database', 'Gene_Symbol'], as_index=False).agg(Pathway_Count=('Pathway_ID', 'nunique'))
print('\nGene coverage depth table:')
print('Shape:', gene_depth.shape)
display(gene_depth.head())
gene_depth_summary = gene_depth.groupby('Database').agg(Mean_Pathways_Per_Gene=('Pathway_Count', 'mean'), Median_Pathways_Per_Gene=('Pathway_Count', 'median'), Min_Pathways_Per_Gene=('Pathway_Count', 'min'), Max_Pathways_Per_Gene=('Pathway_Count', 'max'), SD_Pathways_Per_Gene=('Pathway_Count', 'std')).reset_index()
depth_threshold_summary = gene_depth.assign(In_2plus=lambda x: x['Pathway_Count'] >= 2, In_5plus=lambda x: x['Pathway_Count'] >= 5, In_10plus=lambda x: x['Pathway_Count'] >= 10, In_25plus=lambda x: x['Pathway_Count'] >= 25, In_50plus=lambda x: x['Pathway_Count'] >= 50, In_100plus=lambda x: x['Pathway_Count'] >= 100).groupby('Database').agg(Fraction_Genes_in_2plus_Pathways=('In_2plus', 'mean'), Fraction_Genes_in_5plus_Pathways=('In_5plus', 'mean'), Fraction_Genes_in_10plus_Pathways=('In_10plus', 'mean'), Fraction_Genes_in_25plus_Pathways=('In_25plus', 'mean'), Fraction_Genes_in_50plus_Pathways=('In_50plus', 'mean'), Fraction_Genes_in_100plus_Pathways=('In_100plus', 'mean')).reset_index()
gene_depth_summary = gene_depth_summary.merge(depth_threshold_summary, on='Database', how='left')
top_genes_list = []
for resource in resource_order:
    temp = gene_depth[gene_depth['Database'] == resource].sort_values(['Pathway_Count', 'Gene_Symbol'], ascending=[False, True]).head(20).copy()
    temp['Rank'] = range(1, len(temp) + 1)
    temp = temp[['Database', 'Rank', 'Gene_Symbol', 'Pathway_Count']]
    top_genes_list.append(temp)
top_repeated_genes = pd.concat(top_genes_list, ignore_index=True)
print('\nTop 20 most frequently annotated genes per resource:')
display(top_repeated_genes)
structural_profile = structural_summary.merge(gene_depth_summary, on='Database', how='left')
structural_profile['Database'] = pd.Categorical(structural_profile['Database'], categories=resource_order, ordered=True)
structural_profile = structural_profile.sort_values('Database').reset_index(drop=True)
structural_profile['Pathways_per_1000_Unique_Genes'] = structural_profile['Number_of_Pathways'] / structural_profile['Number_of_Unique_Genes'] * 1000
print('\n' + '=' * 85)
print('STEP 2A — PATHWAY STRUCTURE')
print('=' * 85)
display(structural_profile[['Database', 'Number_of_Pathways', 'Number_of_Unique_Genes', 'Mean_Pathway_Size', 'Median_Pathway_Size', 'Min_Pathway_Size', 'Max_Pathway_Size', 'SD_Pathway_Size']].round(3))
print('\n' + '=' * 85)
print('STEP 2B — GENE COVERAGE DEPTH')
print('=' * 85)
display(structural_profile[['Database', 'Mean_Pathways_Per_Gene', 'Median_Pathways_Per_Gene', 'Min_Pathways_Per_Gene', 'Max_Pathways_Per_Gene', 'Fraction_Genes_in_2plus_Pathways', 'Fraction_Genes_in_5plus_Pathways', 'Fraction_Genes_in_10plus_Pathways', 'Fraction_Genes_in_25plus_Pathways']].round(4))
print('\n' + '=' * 85)
print('STEP 2C — STRUCTURAL INDICATORS')
print('=' * 85)
display(structural_profile[['Database', 'Number_of_Unique_Genes', 'Number_of_Pathways', 'Median_Pathway_Size', 'Fraction_Pathways_10_or_Fewer', 'Fraction_Pathways_25_or_Fewer', 'Fraction_Pathways_Over_100', 'Mean_Pathways_Per_Gene', 'Pathways_per_1000_Unique_Genes']].round(4))
pathway_sizes_file = RESULTS_DIR / 'BRCA_PRISM_step2_pathway_sizes.csv'
structural_profile_file = RESULTS_DIR / 'BRCA_PRISM_step2_structural_profile.csv'
gene_depth_file = RESULTS_DIR / 'BRCA_PRISM_step2_gene_coverage_depth.csv'
gene_depth_summary_file = RESULTS_DIR / 'BRCA_PRISM_step2_gene_coverage_depth_summary.csv'
top_genes_file = RESULTS_DIR / 'BRCA_PRISM_step2_top_repeated_genes.csv'
pathway_sizes.to_csv(pathway_sizes_file, index=False)
structural_profile.to_csv(structural_profile_file, index=False)
gene_depth.to_csv(gene_depth_file, index=False)
gene_depth_summary.to_csv(gene_depth_summary_file, index=False)
top_repeated_genes.to_csv(top_genes_file, index=False)
plot_df = structural_profile.copy()
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(plot_df['Database'].astype(str), plot_df['Number_of_Pathways'])
ax.set_ylabel('Number of pathways / gene sets')
ax.set_xlabel('Resource')
ax.set_title('Number of Pathways Across Resources')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
fig.savefig(FIGURES_DIR / 'BRCA_PRISM_step2_number_of_pathways.png', dpi=300, bbox_inches='tight')
plt.show()
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(plot_df['Database'].astype(str), plot_df['Number_of_Unique_Genes'])
ax.set_ylabel('Number of unique genes')
ax.set_xlabel('Resource')
ax.set_title('Unique Gene Coverage Across Resources')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
fig.savefig(FIGURES_DIR / 'BRCA_PRISM_step2_unique_genes.png', dpi=300, bbox_inches='tight')
plt.show()
box_data = [pathway_sizes.loc[pathway_sizes['Database'] == resource, 'Pathway_Size'].values for resource in resource_order]
fig, ax = plt.subplots(figsize=(9, 6))
ax.boxplot(box_data, tick_labels=resource_order, showfliers=False)
ax.set_yscale('log')
ax.set_ylabel('Genes per pathway (log scale)')
ax.set_xlabel('Resource')
ax.set_title('Pathway-Size Distribution Across Resources')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
fig.savefig(FIGURES_DIR / 'BRCA_PRISM_step2_pathway_size_distribution.png', dpi=300, bbox_inches='tight')
plt.show()
depth_box_data = [gene_depth.loc[gene_depth['Database'] == resource, 'Pathway_Count'].values for resource in resource_order]
fig, ax = plt.subplots(figsize=(9, 6))
ax.boxplot(depth_box_data, tick_labels=resource_order, showfliers=False)
ax.set_yscale('log')
ax.set_ylabel('Number of pathways containing each gene (log scale)')
ax.set_xlabel('Resource')
ax.set_title('Gene Coverage Depth Across Pathway Resources')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
fig.savefig(FIGURES_DIR / 'BRCA_PRISM_step2_gene_coverage_depth.png', dpi=300, bbox_inches='tight')
plt.show()
print('\n' + '=' * 85)
print('FILES SAVED')
print('=' * 85)
print(pathway_sizes_file)
print(structural_profile_file)
print(gene_depth_file)
print(gene_depth_summary_file)
print(top_genes_file)
print('\nFigures saved in:')
print(FIGURES_DIR)
print('\nSTEP 2 structural calculations completed.')
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
FIGURES_DIR = PROJECT_ROOT / 'figures'
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
master = pd.read_csv(MASTER_FILE)
resource_order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB Hallmark']
print('Master shape:', master.shape)
membership = master[['Gene_Symbol', 'Database']].drop_duplicates().assign(Present=1).pivot(index='Gene_Symbol', columns='Database', values='Present').fillna(0).astype(int)
membership = membership[resource_order]
membership['Resource_Count'] = membership.sum(axis=1)
print('Total genes in union:', len(membership))
combination_counts = membership.groupby(resource_order, as_index=False).size().rename(columns={'size': 'Gene_Count'})
combination_counts = combination_counts.sort_values('Gene_Count', ascending=False).reset_index(drop=True)

def readable_combination(row):
    present = [r for r in resource_order if row[r] == 1]
    return ' + '.join(present)
combination_counts['Combination'] = combination_counts.apply(readable_combination, axis=1)
top_n = 15
plot_df = combination_counts.head(top_n).copy().reset_index(drop=True)
plot_df['Rank'] = np.arange(1, len(plot_df) + 1)
rank_key = plot_df[['Rank', 'Combination', 'Gene_Count']].copy()
print('\nTop intersection rank key:')
display(rank_key)
rank_key.to_csv(RESULTS_DIR / 'BRCA_PRISM_gene_overlap_rank_key.csv', index=False)
combination_counts.to_csv(RESULTS_DIR / 'BRCA_PRISM_gene_exact_intersection_counts.csv', index=False)
set_sizes = [membership[r].sum() for r in resource_order]
set_size_df = pd.DataFrame({'Database': resource_order, 'Set_Size': set_sizes})
resource_colors = {'LCPathways': '#1f77b4', 'KEGG': '#ff7f0e', 'Reactome': '#2ca02c', 'WikiPathways': '#d62728', 'MSigDB Hallmark': '#9467bd'}
inactive_color = '#d9d9d9'
bar_color = '#4c78a8'
line_color = '#4d4d4d'
n = len(plot_df)
x = np.arange(n)
y_positions = np.arange(len(resource_order))[::-1]
fig = plt.figure(figsize=(16, 9))
ax_left = fig.add_axes([0.04, 0.12, 0.13, 0.34])
ax_matrix = fig.add_axes([0.27, 0.12, 0.69, 0.34])
ax_bar = fig.add_axes([0.26, 0.54, 0.69, 0.36])
ax_bar.bar(x, plot_df['Gene_Count'], color=bar_color, edgecolor='black', linewidth=0.6)
ax_bar.set_ylabel('Number of genes', fontsize=12)
ax_bar.set_title('Gene overlap across five pathway resources', fontsize=15, pad=10)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels(plot_df['Rank'].astype(str), fontsize=10)
ax_bar.set_xlabel('Intersection rank', fontsize=12)
for i, count in enumerate(plot_df['Gene_Count']):
    ax_bar.text(i, count + max(plot_df['Gene_Count']) * 0.01, str(count), ha='center', va='bottom', fontsize=9)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
for col_idx, (_, row) in enumerate(plot_df.iterrows()):
    active_positions = []
    for resource_idx, resource in enumerate(resource_order):
        y = y_positions[resource_idx]
        if row[resource] == 1:
            ax_matrix.scatter(col_idx, y, s=90, color=resource_colors[resource], edgecolor='black', linewidth=0.4, zorder=3)
            active_positions.append(y)
        else:
            ax_matrix.scatter(col_idx, y, s=55, color=inactive_color, edgecolor='none', zorder=2)
    if len(active_positions) > 1:
        ax_matrix.plot([col_idx, col_idx], [min(active_positions), max(active_positions)], color=line_color, linewidth=1.4, zorder=1)
ax_matrix.set_yticks(y_positions)
ax_matrix.set_yticklabels(resource_order, fontsize=12)
ax_matrix.set_xticks(x)
ax_matrix.set_xticklabels(plot_df['Rank'].astype(str), fontsize=11)
ax_matrix.set_xlabel('Intersection rank', fontsize=12)
ax_matrix.set_xlim(-0.8, n - 0.2)
ax_matrix.set_ylim(-0.6, len(resource_order) - 0.4)
ax_matrix.spines['top'].set_visible(False)
ax_matrix.spines['right'].set_visible(False)
ax_matrix.tick_params(axis='y', pad=2)
for y in y_positions:
    ax_matrix.axhline(y=y, color='#f0f0f0', linewidth=0.8, zorder=0)
ax_left.barh(y_positions, set_sizes, color=bar_color, edgecolor='black', linewidth=0.6)
ax_left.set_yticks([])
ax_left.set_xlabel('Number of unique genes', fontsize=11)
ax_left.set_xlim(max(set_sizes) * 1.05, 0)
for y, count in zip(y_positions, set_sizes):
    ax_left.text(count - 300, y, f'{count:,}', va='center', ha='left', fontsize=10)
ax_left.spines['top'].set_visible(False)
ax_left.spines['right'].set_visible(False)
png_file = FIGURES_DIR / 'BRCA_PRISM_gene_overlap_upset_clean_fullnames.png'
pdf_file = FIGURES_DIR / 'BRCA_PRISM_gene_overlap_upset_clean_fullnames.pdf'
plt.savefig(png_file, dpi=300, bbox_inches='tight')
plt.savefig(pdf_file, bbox_inches='tight')
plt.show()
print('\nSaved figures:')
print(png_file)
print(pdf_file)
support_summary = membership['Resource_Count'].value_counts().sort_index().rename_axis('Number_of_Resources').reset_index(name='Number_of_Genes')
support_summary['Percentage'] = support_summary['Number_of_Genes'] / len(membership) * 100
support_summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_gene_resource_support_summary.csv', index=False)
print('\nGene support summary:')
display(support_summary.round(2))
print('\nINTERPRETATION:')
print('The x-axis labels 1, 2, 3, ..., 15 are INTERSECTION RANKS.')
print('Rank 1 = largest exact gene-overlap combination.')
print('Rank 2 = second largest exact gene-overlap combination.')
print('Rank 3 = third largest exact gene-overlap combination.')
print('Use the saved rank key file to identify each combination:')
print(RESULTS_DIR / 'BRCA_PRISM_gene_overlap_rank_key.csv')
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
master = pd.read_csv(MASTER_FILE)
resource_order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB Hallmark']
print('Master shape:', master.shape)

def analyze_within_resource_overlap(df_resource, resource_name):
    print('\n' + '=' * 80)
    print(resource_name)
    print('=' * 80)
    df = df_resource[['Pathway_ID', 'Pathway_Name', 'Gene_Symbol']].copy()
    df = df.dropna(subset=['Pathway_ID', 'Pathway_Name', 'Gene_Symbol'])
    df['Pathway_ID'] = df['Pathway_ID'].astype(str).str.strip()
    df['Pathway_Name'] = df['Pathway_Name'].astype(str).str.strip()
    df['Gene_Symbol'] = df['Gene_Symbol'].astype(str).str.strip().str.upper()
    invalid_values = {'', 'NAN', 'NONE', 'NA'}
    df = df[~df['Gene_Symbol'].isin(invalid_values)].copy()
    df = df.drop_duplicates(subset=['Pathway_ID', 'Gene_Symbol']).reset_index(drop=True)
    print('Clean associations:', len(df))
    pathway_info = df[['Pathway_ID', 'Pathway_Name']].drop_duplicates(subset=['Pathway_ID']).reset_index(drop=True)
    pathway_ids = pathway_info['Pathway_ID'].tolist()
    genes = sorted(df['Gene_Symbol'].unique().tolist())
    print('Pathways:', len(pathway_ids))
    print('Unique genes:', len(genes))
    pathway_index = {pathway_id: i for i, pathway_id in enumerate(pathway_ids)}
    gene_index = {gene: i for i, gene in enumerate(genes)}
    rows = df['Pathway_ID'].map(pathway_index).to_numpy()
    cols = df['Gene_Symbol'].map(gene_index).to_numpy()
    data = np.ones(len(df), dtype=np.int32)
    X = csr_matrix((data, (rows, cols)), shape=(len(pathway_ids), len(genes)), dtype=np.int32)
    X.data[:] = 1
    pathway_sizes = np.asarray(X.sum(axis=1)).ravel()
    intersections = X @ X.T
    intersections.setdiag(0)
    intersections.eliminate_zeros()
    best_rows = []
    for i, pathway_id in enumerate(pathway_ids):
        start = intersections.indptr[i]
        end = intersections.indptr[i + 1]
        candidate_indices = intersections.indices[start:end]
        shared_counts = intersections.data[start:end]
        best_jaccard = 0.0
        best_idx = None
        best_shared = 0
        for j, shared in zip(candidate_indices, shared_counts):
            union = pathway_sizes[i] + pathway_sizes[j] - shared
            if union <= 0:
                continue
            jaccard = float(shared) / float(union)
            if jaccard > best_jaccard:
                best_jaccard = jaccard
                best_idx = j
                best_shared = shared
        if best_idx is None:
            best_rows.append({'Database': resource_name, 'Pathway_ID': pathway_id, 'Pathway_Name': pathway_info.loc[i, 'Pathway_Name'], 'Pathway_Size': int(pathway_sizes[i]), 'Best_Match_ID': None, 'Best_Match_Name': None, 'Best_Match_Size': np.nan, 'Shared_Genes': 0, 'Max_Jaccard': 0.0})
        else:
            best_rows.append({'Database': resource_name, 'Pathway_ID': pathway_id, 'Pathway_Name': pathway_info.loc[i, 'Pathway_Name'], 'Pathway_Size': int(pathway_sizes[i]), 'Best_Match_ID': pathway_ids[best_idx], 'Best_Match_Name': pathway_info.loc[best_idx, 'Pathway_Name'], 'Best_Match_Size': int(pathway_sizes[best_idx]), 'Shared_Genes': int(best_shared), 'Max_Jaccard': float(best_jaccard)})
    best_df = pd.DataFrame(best_rows)
    summary = {'Database': resource_name, 'Number_of_Pathways': len(best_df), 'Mean_Max_Jaccard': best_df['Max_Jaccard'].mean(), 'Median_Max_Jaccard': best_df['Max_Jaccard'].median(), 'Max_of_Max_Jaccard': best_df['Max_Jaccard'].max(), 'Fraction_Pathways_MaxJaccard_ge_0.25': (best_df['Max_Jaccard'] >= 0.25).mean(), 'Fraction_Pathways_MaxJaccard_ge_0.50': (best_df['Max_Jaccard'] >= 0.5).mean(), 'Fraction_Pathways_MaxJaccard_ge_0.75': (best_df['Max_Jaccard'] >= 0.75).mean(), 'Fraction_Pathways_MaxJaccard_ge_0.90': (best_df['Max_Jaccard'] >= 0.9).mean()}
    print('Mean maximum Jaccard:', round(summary['Mean_Max_Jaccard'], 4))
    print('Median maximum Jaccard:', round(summary['Median_Max_Jaccard'], 4))
    print('Fraction >= 0.50:', round(summary['Fraction_Pathways_MaxJaccard_ge_0.50'], 4))
    return (best_df, summary)
all_best_matches = []
all_summaries = []
for resource in resource_order:
    resource_df = master[master['Database'] == resource].copy()
    best_df, summary = analyze_within_resource_overlap(resource_df, resource)
    all_best_matches.append(best_df)
    all_summaries.append(summary)
within_overlap = pd.concat(all_best_matches, ignore_index=True)
redundancy_summary = pd.DataFrame(all_summaries)
redundancy_summary['Database'] = pd.Categorical(redundancy_summary['Database'], categories=resource_order, ordered=True)
redundancy_summary = redundancy_summary.sort_values('Database').reset_index(drop=True)
print('\n' + '=' * 95)
print('WITHIN-RESOURCE PATHWAY OVERLAP SUMMARY')
print('=' * 95)
display(redundancy_summary.round(4))
top_pairs_list = []
for resource in resource_order:
    tmp = within_overlap[within_overlap['Database'] == resource].sort_values('Max_Jaccard', ascending=False).head(10).copy()
    top_pairs_list.append(tmp)
top_pairs = pd.concat(top_pairs_list, ignore_index=True)
print('\nTop 10 strongest within-resource pathway overlaps:')
display(top_pairs[['Database', 'Pathway_Name', 'Best_Match_Name', 'Pathway_Size', 'Best_Match_Size', 'Shared_Genes', 'Max_Jaccard']].round(4))
within_overlap_file = RESULTS_DIR / 'BRCA_PRISM_step2_within_resource_pathway_overlap.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_step2_within_resource_overlap_summary.csv'
top_pairs_file = RESULTS_DIR / 'BRCA_PRISM_step2_top_within_resource_overlaps.csv'
within_overlap.to_csv(within_overlap_file, index=False)
redundancy_summary.to_csv(summary_file, index=False)
top_pairs.to_csv(top_pairs_file, index=False)
plot_data = [within_overlap.loc[within_overlap['Database'] == resource, 'Max_Jaccard'].values for resource in resource_order]
fig, ax = plt.subplots(figsize=(9, 6))
ax.boxplot(plot_data, tick_labels=resource_order, showfliers=False)
ax.set_ylabel('Maximum within-resource Jaccard similarity')
ax.set_xlabel('Pathway resource')
ax.set_title('Within-Resource Pathway Gene-Content Overlap')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
figure_file = FIGURES_DIR / 'BRCA_PRISM_step2_within_resource_pathway_overlap.png'
plt.savefig(figure_file, dpi=300, bbox_inches='tight')
plt.show()
print('\nSaved:')
print(within_overlap_file)
print(summary_file)
print(top_pairs_file)
print(figure_file)
print('\nSTEP 2 within-resource pathway overlap analysis completed.')
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
PROJECT_ROOT = Path(__file__).resolve().parents[1]
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
RESULTS_DIR = PROJECT_ROOT / 'results'
FIGURES_DIR = PROJECT_ROOT / 'figures'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
master = pd.read_csv(MASTER_FILE)
resource_order = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB Hallmark']
master = master.dropna(subset=['Database', 'Gene_Symbol']).copy()
master['Database'] = master['Database'].astype(str).str.strip()
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master = master[master['Database'].isin(resource_order)].copy()
print('Master shape:', master.shape)
print('Resources:', master['Database'].unique().tolist())
gene_sets = {}
for resource in resource_order:
    genes = set(master.loc[master['Database'] == resource, 'Gene_Symbol'].dropna().unique())
    gene_sets[resource] = genes
    print(f'{resource:<18}: {len(genes):,} unique genes')
pairwise_rows = []
for i, resource_a in enumerate(resource_order):
    for j in range(i + 1, len(resource_order)):
        resource_b = resource_order[j]
        A = gene_sets[resource_a]
        B = gene_sets[resource_b]
        intersection = A & B
        union = A | B
        shared = len(intersection)
        jaccard = shared / len(union) if len(union) > 0 else np.nan
        pairwise_rows.append({'Resource_A': resource_a, 'Resource_B': resource_b, 'Genes_A': len(A), 'Genes_B': len(B), 'Shared_Genes': shared, 'Union_Genes': len(union), 'Jaccard': jaccard})
pairwise_df = pd.DataFrame(pairwise_rows)
print('\n' + '=' * 90)
print('PAIRWISE SHARED GENES AND JACCARD SIMILARITY')
print('=' * 90)
display(pairwise_df.round(4))
jaccard_matrix = pd.DataFrame(index=resource_order, columns=resource_order, dtype=float)
shared_matrix = pd.DataFrame(index=resource_order, columns=resource_order, dtype=int)
for resource_a in resource_order:
    for resource_b in resource_order:
        A = gene_sets[resource_a]
        B = gene_sets[resource_b]
        shared = len(A & B)
        union = len(A | B)
        shared_matrix.loc[resource_a, resource_b] = shared
        jaccard_matrix.loc[resource_a, resource_b] = shared / union if union > 0 else np.nan
print('\nPairwise shared-gene matrix:')
display(shared_matrix)
print('\nPairwise Jaccard matrix:')
display(jaccard_matrix.round(4))
directional_rows = []
directional_matrix = pd.DataFrame(index=resource_order, columns=resource_order, dtype=float)
for source in resource_order:
    for target in resource_order:
        source_genes = gene_sets[source]
        target_genes = gene_sets[target]
        shared = len(source_genes & target_genes)
        coverage = shared / len(source_genes) if len(source_genes) > 0 else np.nan
        directional_matrix.loc[source, target] = coverage
        if source != target:
            directional_rows.append({'Source_Resource': source, 'Target_Resource': target, 'Source_Gene_Count': len(source_genes), 'Shared_Genes': shared, 'Directional_Coverage': coverage})
directional_df = pd.DataFrame(directional_rows)
print('\n' + '=' * 90)
print('DIRECTIONAL GENE COVERAGE')
print('=' * 90)
display(directional_df.round(4))
print('\nDirectional coverage matrix:')
display(directional_matrix.round(4))
membership = master[['Gene_Symbol', 'Database']].drop_duplicates().assign(Present=1).pivot(index='Gene_Symbol', columns='Database', values='Present').fillna(0).astype(int)
membership = membership[resource_order]
membership['Resource_Count'] = membership[resource_order].sum(axis=1)
print('\nTotal unique genes in union:')
print(f'{len(membership):,}')
support_summary = membership['Resource_Count'].value_counts().sort_index().rename_axis('Number_of_Resources').reset_index(name='Number_of_Genes')
support_summary['Percentage_of_Union'] = support_summary['Number_of_Genes'] / len(membership) * 100
print('\n' + '=' * 90)
print('GENE SUPPORT ACROSS 1–5 RESOURCES')
print('=' * 90)
display(support_summary.round(2))
consensus_genes = membership[membership['Resource_Count'] >= 3].copy()
consensus_gene_list = consensus_genes.reset_index()
core_genes = membership[membership['Resource_Count'] == 5].copy()
core_gene_list = core_genes.reset_index()
print('\n' + '=' * 90)
print('GLOBAL CONSENSUS AND CORE GENE SETS')
print('=' * 90)
print('Consensus genes (>=3 resources):', f'{len(consensus_gene_list):,}')
print('Core genes (all 5 resources):', f'{len(core_gene_list):,}')
print('Consensus genes as % of union:', round(len(consensus_gene_list) / len(membership) * 100, 2), '%')
print('Core genes as % of union:', round(len(core_gene_list) / len(membership) * 100, 2), '%')
resource_specific_rows = []
resource_specific_lists = []
for resource in resource_order:
    other_resources = [r for r in resource_order if r != resource]
    specific_mask = (membership[resource] == 1) & (membership[other_resources].sum(axis=1) == 0)
    specific_genes = membership[specific_mask].reset_index()[['Gene_Symbol']].copy()
    specific_genes['Database'] = resource
    specific_genes = specific_genes[['Database', 'Gene_Symbol']]
    resource_specific_lists.append(specific_genes)
    resource_specific_rows.append({'Database': resource, 'Resource_Specific_Genes': len(specific_genes), 'Total_Resource_Genes': len(gene_sets[resource]), 'Fraction_of_Resource_Specific': len(specific_genes) / len(gene_sets[resource])})
resource_specific_summary = pd.DataFrame(resource_specific_rows)
resource_specific_genes = pd.concat(resource_specific_lists, ignore_index=True)
print('\n' + '=' * 90)
print('RESOURCE-SPECIFIC GENES')
print('=' * 90)
display(resource_specific_summary.round(4))
pairwise_df.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_pairwise_gene_overlap.csv', index=False)
shared_matrix.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_shared_gene_matrix.csv')
jaccard_matrix.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_gene_jaccard_matrix.csv')
directional_df.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_directional_gene_coverage.csv', index=False)
directional_matrix.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_directional_gene_coverage_matrix.csv')
membership.reset_index().to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_gene_resource_membership.csv', index=False)
support_summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_gene_support_summary.csv', index=False)
consensus_gene_list.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_global_consensus_genes_ge3.csv', index=False)
core_gene_list.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_core_genes_all5.csv', index=False)
resource_specific_summary.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_resource_specific_gene_summary.csv', index=False)
resource_specific_genes.to_csv(RESULTS_DIR / 'BRCA_PRISM_step3_resource_specific_genes.csv', index=False)
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(jaccard_matrix.values, aspect='auto')
ax.set_xticks(np.arange(len(resource_order)))
ax.set_yticks(np.arange(len(resource_order)))
ax.set_xticklabels(resource_order, rotation=35, ha='right')
ax.set_yticklabels(resource_order)
for i in range(len(resource_order)):
    for j in range(len(resource_order)):
        ax.text(j, i, f'{jaccard_matrix.iloc[i, j]:.2f}', ha='center', va='center', fontsize=9)
ax.set_title('Pairwise Gene-Level Jaccard Similarity')
fig.colorbar(im, ax=ax, label='Jaccard similarity')
plt.tight_layout()
jaccard_fig = FIGURES_DIR / 'BRCA_PRISM_step3_gene_jaccard_heatmap.png'
plt.savefig(jaccard_fig, dpi=300, bbox_inches='tight')
plt.show()
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(directional_matrix.values, aspect='auto')
ax.set_xticks(np.arange(len(resource_order)))
ax.set_yticks(np.arange(len(resource_order)))
ax.set_xticklabels(resource_order, rotation=35, ha='right')
ax.set_yticklabels(resource_order)
for i in range(len(resource_order)):
    for j in range(len(resource_order)):
        ax.text(j, i, f'{directional_matrix.iloc[i, j]:.2f}', ha='center', va='center', fontsize=9)
ax.set_xlabel('Target resource')
ax.set_ylabel('Source resource')
ax.set_title('Directional Gene Coverage Across Pathway Resources')
fig.colorbar(im, ax=ax, label='Fraction of source genes covered')
plt.tight_layout()
directional_fig = FIGURES_DIR / 'BRCA_PRISM_step3_directional_gene_coverage_heatmap.png'
plt.savefig(directional_fig, dpi=300, bbox_inches='tight')
plt.show()
fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(support_summary['Number_of_Resources'].astype(str), support_summary['Number_of_Genes'])
ax.set_xlabel('Number of pathway resources containing a gene')
ax.set_ylabel('Number of genes')
ax.set_title('Cross-Resource Support for Genes')
for i, row in support_summary.iterrows():
    ax.text(i, row['Number_of_Genes'], f"{int(row['Number_of_Genes']):,}", ha='center', va='bottom')
plt.tight_layout()
support_fig = FIGURES_DIR / 'BRCA_PRISM_step3_gene_support_distribution.png'
plt.savefig(support_fig, dpi=300, bbox_inches='tight')
plt.show()
print('\n' + '=' * 90)
print('STEP 3 COMPLETE')
print('=' * 90)
print('\nSaved major figures:')
print(jaccard_fig)
print(directional_fig)
print(support_fig)
print('\nMain outputs to send me:')
print('1. Pairwise shared genes and Jaccard table')
print('2. Directional coverage matrix')
print('3. Gene support summary')
print('4. Resource-specific gene summary')
print('5. Consensus and core gene counts')
