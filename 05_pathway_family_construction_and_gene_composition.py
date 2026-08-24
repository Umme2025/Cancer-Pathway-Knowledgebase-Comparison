from pathlib import Path
import pandas as pd
import networkx as nx
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
STEP5_FILE = RESULTS_DIR / 'BRCA_PRISM_step5_refined_provisional_pathway_relationships.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
relations = pd.read_csv(STEP5_FILE)
master = pd.read_csv(MASTER_FILE)
print('Step 5 relationships:', relations.shape)
print('Master:', master.shape)
equiv = relations[relations['Provisional_Relationship'] == 'Equivalent candidate'].copy()
print('Directional equivalent mappings:', len(equiv))
G = nx.Graph()
for _, row in equiv.iterrows():
    source_node = (str(row['Source_Resource']), str(row['Source_Pathway_ID']))
    target_node = (str(row['Target_Resource']), str(row['Target_Pathway_ID']))
    G.add_edge(source_node, target_node)
components = list(nx.connected_components(G))
components = sorted(components, key=len, reverse=True)
print('Initial BRCA-PRISM families:', len(components))
family_rows = []
for i, component in enumerate(components, start=1):
    family_id = f'BRCA_PRISM_F{i:04d}'
    for database, pathway_id in sorted(component):
        subset = master[(master['Database'] == database) & (master['Pathway_ID'].astype(str) == str(pathway_id))]
        if len(subset) > 0:
            pathway_name = subset['Pathway_Name'].iloc[0]
            gene_count = subset['Gene_Symbol'].nunique()
        else:
            pathway_name = ''
            gene_count = 0
        family_rows.append({'Family_ID': family_id, 'Family_Name': '', 'Database': database, 'Pathway_ID': pathway_id, 'Pathway_Name': pathway_name, 'Gene_Count': gene_count, 'Family_Membership_Evidence': 'Equivalent candidate', 'Original_Pathway_Preserved': True})
family_layer = pd.DataFrame(family_rows)
family_summary = family_layer.groupby('Family_ID').agg(Number_of_Pathways=('Pathway_ID', 'count'), Number_of_Resources=('Database', 'nunique'), Resources=('Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
print('\nInitial pathway-family mapping layer:')
display(family_layer.head(50))
print('\nFamily summary:')
display(family_summary.head(30))
family_file = RESULTS_DIR / 'BRCA_PRISM_step6A_initial_pathway_family_layer.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_step6A_initial_family_summary.csv'
family_layer.to_csv(family_file, index=False)
family_summary.to_csv(summary_file, index=False)
print('\nSTEP 6A COMPLETE')
print('\nSaved:')
print(family_file)
print(summary_file)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6A_initial_pathway_family_layer.csv'
RELATION_FILE = RESULTS_DIR / 'BRCA_PRISM_step5_refined_provisional_pathway_relationships.csv'
families = pd.read_csv(FAMILY_FILE)
relations = pd.read_csv(RELATION_FILE)
print('Initial family pathways:', len(families))
print('Step 5 relationships:', len(relations))
family_lookup = {}
for _, row in families.iterrows():
    key = (str(row['Database']), str(row['Pathway_ID']))
    family_lookup[key] = row['Family_ID']
candidate_relations = relations[relations['Provisional_Relationship'].isin(['Hierarchical / isPartOf candidate', 'Related candidate'])].copy()
rows = []
for _, row in candidate_relations.iterrows():
    source_key = (str(row['Source_Resource']), str(row['Source_Pathway_ID']))
    target_key = (str(row['Target_Resource']), str(row['Target_Pathway_ID']))
    source_family = family_lookup.get(source_key)
    target_family = family_lookup.get(target_key)
    if source_family is not None and target_family is None:
        rows.append({'Family_ID': source_family, 'Anchor_Database': row['Source_Resource'], 'Anchor_Pathway_ID': row['Source_Pathway_ID'], 'Anchor_Pathway_Name': row['Source_Pathway_Name'], 'Candidate_Database': row['Target_Resource'], 'Candidate_Pathway_ID': row['Target_Pathway_ID'], 'Candidate_Pathway_Name': row['Target_Pathway_Name'], 'Relationship': row['Provisional_Relationship'], 'Shared_Genes': row['Shared_Genes'], 'Jaccard': row['Jaccard'], 'Overlap_Coefficient': row['Overlap_Coefficient'], 'Name_Similarity': row['Name_Similarity'], 'Reciprocal_Best_Match': row['Reciprocal_Best_Match'], 'Containment_Direction': row['Containment_Direction']})
    elif target_family is not None and source_family is None:
        rows.append({'Family_ID': target_family, 'Anchor_Database': row['Target_Resource'], 'Anchor_Pathway_ID': row['Target_Pathway_ID'], 'Anchor_Pathway_Name': row['Target_Pathway_Name'], 'Candidate_Database': row['Source_Resource'], 'Candidate_Pathway_ID': row['Source_Pathway_ID'], 'Candidate_Pathway_Name': row['Source_Pathway_Name'], 'Relationship': row['Provisional_Relationship'], 'Shared_Genes': row['Shared_Genes'], 'Jaccard': row['Jaccard'], 'Overlap_Coefficient': row['Overlap_Coefficient'], 'Name_Similarity': row['Name_Similarity'], 'Reciprocal_Best_Match': row['Reciprocal_Best_Match'], 'Containment_Direction': row['Containment_Direction']})
expansion_candidates = pd.DataFrame(rows)
expansion_candidates = expansion_candidates.drop_duplicates(subset=['Family_ID', 'Candidate_Database', 'Candidate_Pathway_ID'])
expansion_candidates = expansion_candidates.sort_values(['Family_ID', 'Relationship', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity'], ascending=[True, True, False, False, False]).reset_index(drop=True)
expansion_candidates['Add_to_Family'] = ''
expansion_candidates['Final_Relationship'] = ''
expansion_candidates['Review_Notes'] = ''
print('\n' + '=' * 90)
print('STEP 6B FAMILY-EXPANSION CANDIDATES')
print('=' * 90)
print('Total candidate pathways:', len(expansion_candidates))
print('Families with expansion candidates:', expansion_candidates['Family_ID'].nunique())
display(expansion_candidates.head(50))
family_candidate_summary = expansion_candidates.groupby('Family_ID').agg(Number_of_Candidate_Pathways=('Candidate_Pathway_ID', 'nunique'), Candidate_Resources=('Candidate_Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
print('\nCandidate summary by family:')
display(family_candidate_summary.head(31))
OUTPUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step6B_family_expansion_review.csv'
SUMMARY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6B_family_expansion_summary.csv'
expansion_candidates.to_csv(OUTPUT_FILE, index=False)
family_candidate_summary.to_csv(SUMMARY_FILE, index=False)
print('\nSTEP 6B COMPLETE')
print('\nSaved:')
print(OUTPUT_FILE)
print(SUMMARY_FILE)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
INPUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step6B_family_expansion_review.csv'
df = pd.read_csv(INPUT_FILE)
selected_families = ['BRCA_PRISM_F0002', 'BRCA_PRISM_F0006', 'BRCA_PRISM_F0008', 'BRCA_PRISM_F0009', 'BRCA_PRISM_F0010', 'BRCA_PRISM_F0011', 'BRCA_PRISM_F0012', 'BRCA_PRISM_F0013', 'BRCA_PRISM_F0015', 'BRCA_PRISM_F0016', 'BRCA_PRISM_F0018']
review = df[df['Family_ID'].isin(selected_families)].copy()
review = review.sort_values(['Family_ID', 'Relationship', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity'], ascending=[True, True, False, False, False]).reset_index(drop=True)
review['Add_to_Family'] = ''
review['Final_Relationship'] = ''
review['Review_Notes'] = ''
print('BRCA-focused candidates:', len(review))
print('Families represented:', review['Family_ID'].nunique())
display(review[['Family_ID', 'Anchor_Database', 'Anchor_Pathway_Name', 'Candidate_Database', 'Candidate_Pathway_ID', 'Candidate_Pathway_Name', 'Relationship', 'Shared_Genes', 'Jaccard', 'Overlap_Coefficient', 'Name_Similarity', 'Reciprocal_Best_Match', 'Add_to_Family', 'Final_Relationship', 'Review_Notes']])
OUTPUT_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_BRCA_family_manual_review.csv'
review.to_csv(OUTPUT_FILE, index=False)
print('\nSaved:')
print(OUTPUT_FILE)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
INITIAL_FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6A_initial_pathway_family_layer.csv'
REVIEW_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_BRCA_family_manual_review_COMPLETED.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
families = pd.read_csv(INITIAL_FAMILY_FILE)
review = pd.read_csv(REVIEW_FILE)
master = pd.read_csv(MASTER_FILE)
approved = review[review['Add_to_Family'].astype(str).str.strip().str.lower() == 'yes'].copy()
print('Initial family rows:', len(families))
print('Approved expansion pathways:', len(approved))
new_rows = []
for _, row in approved.iterrows():
    db = str(row['Candidate_Database'])
    pid = str(row['Candidate_Pathway_ID'])
    subset = master[(master['Database'] == db) & (master['Pathway_ID'].astype(str) == pid)]
    if len(subset) == 0:
        print('WARNING not found:', db, pid)
        continue
    pathway_name = subset['Pathway_Name'].iloc[0]
    gene_count = subset['Gene_Symbol'].nunique()
    new_rows.append({'Family_ID': row['Family_ID'], 'Family_Name': '', 'Database': db, 'Pathway_ID': pid, 'Pathway_Name': pathway_name, 'Gene_Count': gene_count, 'Family_Membership_Evidence': row['Final_Relationship'], 'Original_Pathway_Preserved': True})
approved_layer = pd.DataFrame(new_rows)
final_family_layer = pd.concat([families, approved_layer], ignore_index=True)
final_family_layer = final_family_layer.drop_duplicates(subset=['Family_ID', 'Database', 'Pathway_ID']).reset_index(drop=True)
final_summary = final_family_layer.groupby('Family_ID').agg(Number_of_Pathways=('Pathway_ID', 'count'), Number_of_Resources=('Database', 'nunique'), Resources=('Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
print('\nFinal family-layer rows:', len(final_family_layer))
print('Number of families:', final_family_layer['Family_ID'].nunique())
display(final_family_layer.head(60))
print('\nFinal family summary:')
display(final_summary.head(31))
FINAL_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
SUMMARY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_family_summary.csv'
final_family_layer.to_csv(FINAL_FILE, index=False)
final_summary.to_csv(SUMMARY_FILE, index=False)
print('\nSTEP 6 COMPLETE')
print('\nSaved:')
print(FINAL_FILE)
print(SUMMARY_FILE)
print('\nSend me:')
print('1. Final family-layer rows')
print('2. Number of families')
print('3. First 20 rows of final_summary')
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
FAMILY_FILE = RESULTS_DIR / 'BRCA_PRISM_step6C_final_reviewed_pathway_family_layer.csv'
MASTER_FILE = PROJECT_ROOT / 'data_processed' / 'BRCA_PRISM_all_pathway_resources_master.csv'
family_layer = pd.read_csv(FAMILY_FILE)
master = pd.read_csv(MASTER_FILE)
print('Family layer:', family_layer.shape)
print('Master:', master.shape)
master = master.dropna(subset=['Database', 'Pathway_ID', 'Gene_Symbol']).copy()
master['Gene_Symbol'] = master['Gene_Symbol'].astype(str).str.strip().str.upper()
master['Pathway_ID'] = master['Pathway_ID'].astype(str)
family_layer['Pathway_ID'] = family_layer['Pathway_ID'].astype(str)
family_genes = family_layer.merge(master[['Database', 'Pathway_ID', 'Gene_Symbol']], on=['Database', 'Pathway_ID'], how='left')
family_genes = family_genes.dropna(subset=['Gene_Symbol'])
family_genes = family_genes.drop_duplicates(subset=['Family_ID', 'Database', 'Pathway_ID', 'Gene_Symbol'])
print('Family-gene associations:', len(family_genes))
gene_support = family_genes[['Family_ID', 'Database', 'Gene_Symbol']].drop_duplicates().groupby(['Family_ID', 'Gene_Symbol']).agg(Resource_Support=('Database', 'nunique'), Resources=('Database', lambda x: '; '.join(sorted(set(x))))).reset_index()
family_resource_count = family_layer.groupby('Family_ID')['Database'].nunique().to_dict()

def classify_gene(row):
    family_id = row['Family_ID']
    support = row['Resource_Support']
    n_resources = family_resource_count[family_id]
    if support == n_resources:
        return 'Core'
    elif support >= 2:
        return 'Consensus'
    else:
        return 'Resource-specific'
gene_support['Gene_Category'] = gene_support.apply(classify_gene, axis=1)
family_summary_rows = []
for family_id, group in gene_support.groupby('Family_ID'):
    union_n = group['Gene_Symbol'].nunique()
    core_n = (group['Gene_Category'] == 'Core').sum()
    consensus_n = (group['Resource_Support'] >= 2).sum()
    specific_n = (group['Gene_Category'] == 'Resource-specific').sum()
    family_summary_rows.append({'Family_ID': family_id, 'Number_of_Resources': family_resource_count[family_id], 'Union_Genes': union_n, 'Consensus_Genes_ge2': consensus_n, 'Core_Genes_all_resources': core_n, 'Resource_Specific_Genes': specific_n})
family_gene_summary = pd.DataFrame(family_summary_rows)
resource_specific = gene_support[gene_support['Gene_Category'] == 'Resource-specific'].copy()
resource_specific_summary = resource_specific.groupby(['Family_ID', 'Resources']).size().reset_index(name='Number_of_Resource_Specific_Genes')
print('\n' + '=' * 100)
print('STEP 7 FAMILY GENE SUMMARY')
print('=' * 100)
display(family_gene_summary.sort_values('Family_ID'))
print('\nExample gene-support rows:')
display(gene_support.head(50))
SUPPORT_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
SUMMARY_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_summary.csv'
SPECIFIC_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_resource_specific_gene_summary.csv'
gene_support.to_csv(SUPPORT_FILE, index=False)
family_gene_summary.to_csv(SUMMARY_FILE, index=False)
resource_specific_summary.to_csv(SPECIFIC_FILE, index=False)
print('\nSTEP 7 COMPLETE')
print('\nSaved:')
print(SUPPORT_FILE)
print(SUMMARY_FILE)
print(SPECIFIC_FILE)
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
TCGA_FILE = PROJECT_ROOT / 'data_raw' / 'TCGA_BRCA' / 'tcga_brca_as_validation.csv'
GENE_SUPPORT_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_support.csv'
SUMMARY_FILE = RESULTS_DIR / 'BRCA_PRISM_step7_family_gene_summary.csv'
tcga = pd.read_csv(TCGA_FILE, nrows=2)
gene_support = pd.read_csv(GENE_SUPPORT_FILE)
family_summary = pd.read_csv(SUMMARY_FILE)
print('TCGA columns:', len(tcga.columns))
print('Step 7 gene-support rows:', len(gene_support))
rna_cols = [c for c in tcga.columns if str(c).lower().endswith('_rna')]
cnv_cols = [c for c in tcga.columns if str(c).lower().endswith('_cnv')]
rna_genes = {c[:-4].strip().upper() for c in rna_cols}
cnv_genes = {c[:-4].strip().upper() for c in cnv_cols}
tcga_union_genes = rna_genes | cnv_genes
tcga_both_genes = rna_genes & cnv_genes
print('\nTCGA RNA genes:', len(rna_genes))
print('TCGA CNV genes:', len(cnv_genes))
print('TCGA RNA∪CNV measured genes:', len(tcga_union_genes))
print('Genes measured in BOTH RNA and CNV:', len(tcga_both_genes))
gene_support['Gene_Symbol'] = gene_support['Gene_Symbol'].astype(str).str.strip().str.upper()
gene_support['Measured_in_TCGA_RNA'] = gene_support['Gene_Symbol'].isin(rna_genes)
gene_support['Measured_in_TCGA_CNV'] = gene_support['Gene_Symbol'].isin(cnv_genes)
gene_support['Measured_in_TCGA_RNA_or_CNV'] = gene_support['Gene_Symbol'].isin(tcga_union_genes)
gene_support['Measured_in_Both_RNA_and_CNV'] = gene_support['Gene_Symbol'].isin(tcga_both_genes)
rows = []
for family_id, group in gene_support.groupby('Family_ID'):
    union_total = len(group)
    measured_union = group['Measured_in_TCGA_RNA_or_CNV'].sum()
    both_omics = group['Measured_in_Both_RNA_and_CNV'].sum()
    consensus_group = group[group['Resource_Support'] >= 2]
    core_group = group[group['Gene_Category'] == 'Core']
    consensus_total = len(consensus_group)
    consensus_measured = consensus_group['Measured_in_TCGA_RNA_or_CNV'].sum()
    core_total = len(core_group)
    core_measured = core_group['Measured_in_TCGA_RNA_or_CNV'].sum()
    rows.append({'Family_ID': family_id, 'Union_Genes': union_total, 'TCGA_Measured_Union_Genes': int(measured_union), 'Union_Retention': measured_union / union_total if union_total > 0 else 0, 'Genes_Measured_in_Both_RNA_CNV': int(both_omics), 'Consensus_Genes_ge2': consensus_total, 'TCGA_Measured_Consensus_Genes': int(consensus_measured), 'Consensus_Retention': consensus_measured / consensus_total if consensus_total > 0 else 0, 'Core_Genes': core_total, 'TCGA_Measured_Core_Genes': int(core_measured), 'Core_Retention': core_measured / core_total if core_total > 0 else 0})
tcga_family_coverage = pd.DataFrame(rows)
tcga_family_coverage = tcga_family_coverage.merge(family_summary[['Family_ID', 'Number_of_Resources']], on='Family_ID', how='left')
print('\n' + '=' * 100)
print('STEP 8 TCGA-BRCA FAMILY COVERAGE')
print('=' * 100)
display(tcga_family_coverage.sort_values('Family_ID').round(4))
print('\nOverall mean family union retention:', round(tcga_family_coverage['Union_Retention'].mean(), 4))
print('Overall median family union retention:', round(tcga_family_coverage['Union_Retention'].median(), 4))
print('Mean consensus retention:', round(tcga_family_coverage['Consensus_Retention'].mean(), 4))
print('Mean core retention:', round(tcga_family_coverage['Core_Retention'].mean(), 4))
GENE_LEVEL_FILE = RESULTS_DIR / 'BRCA_PRISM_step8_TCGA_family_gene_measurement_status.csv'
COVERAGE_FILE = RESULTS_DIR / 'BRCA_PRISM_step8_TCGA_family_coverage.csv'
gene_support.to_csv(GENE_LEVEL_FILE, index=False)
tcga_family_coverage.to_csv(COVERAGE_FILE, index=False)
print('\nSTEP 8 COMPLETE')
print('\nSaved:')
print(GENE_LEVEL_FILE)
print(COVERAGE_FILE)
print('\nSend me:')
print('1. TCGA RNA genes / CNV genes / union genes')
print('2. Full STEP 8 TCGA-BRCA FAMILY COVERAGE table')
print('3. Mean/median union retention')
print('4. Mean consensus retention and mean core retention')
