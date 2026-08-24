from pathlib import Path
BASE_DIR = Path.cwd()
RAW_DIR = BASE_DIR / 'data_raw'
print('Project folder:', BASE_DIR)
print('Raw data folder:', RAW_DIR)
print()
for folder in RAW_DIR.iterdir():
    if folder.is_dir():
        print('=' * 60)
        print('DATABASE:', folder.name)
        print('=' * 60)
        files = list(folder.iterdir())
        print('Number of files:', len(files))
        for f in files[:20]:
            print(' -', f.name)
        if len(files) > 20:
            print(f' ... and {len(files) - 20} more files')
        print()
import json
import pandas as pd
from pathlib import Path
RAW_DIR = Path.cwd() / 'data_raw'
hallmark_dir = RAW_DIR / 'Hallmarks_Pathways'
hallmark_file = sorted(hallmark_dir.glob('*.json'))[0]
print('\n' + '=' * 70)
print('HALLMARK')
print('=' * 70)
print('Example file:', hallmark_file.name)
with open(hallmark_file, 'r', encoding='utf-8') as f:
    hallmark_obj = json.load(f)
print('Top-level keys:', list(hallmark_obj.keys()))
first_pathway = next(iter(hallmark_obj))
print('Example pathway:', first_pathway)
print('Fields:', hallmark_obj[first_pathway].keys())
print('Number of genes:', len(hallmark_obj[first_pathway]['geneSymbols']))
print('First 10 genes:', hallmark_obj[first_pathway]['geneSymbols'][:10])
kegg_file = RAW_DIR / 'KEGG_Homo_sapiens_final_database' / 'KEGG_Homo_sapiens_pathway_gene_database.csv'
print('\n' + '=' * 70)
print('KEGG')
print('=' * 70)
kegg = pd.read_csv(kegg_file)
print('Shape:', kegg.shape)
print('Columns:', list(kegg.columns))
print(kegg.head())
lc_file = RAW_DIR / 'LCPathways' / '41568_2020_240_MOESM4_ESM.csv'
print('\n' + '=' * 70)
print('LCPATHWAYS')
print('=' * 70)
lc = pd.read_csv(lc_file)
print('Shape:', lc.shape)
print('Columns:', list(lc.columns))
print(lc.head())
reactome_file = RAW_DIR / 'ReactomePathways.gmt' / 'ReactomePathways.gmt'
print('\n' + '=' * 70)
print('REACTOME')
print('=' * 70)
with open(reactome_file, 'r', encoding='utf-8') as f:
    first_line = f.readline().rstrip('\n')
parts = first_line.split('\t')
print('Number of fields in first row:', len(parts))
print('Pathway field 1:', parts[0])
print('Pathway field 2:', parts[1] if len(parts) > 1 else None)
print('First 10 genes:', parts[2:12])
wiki_file = RAW_DIR / 'wikipathways' / 'wikipathways-20230210-gmt-Homo_sapiens.gmt'
print('\n' + '=' * 70)
print('WIKIPATHWAYS')
print('=' * 70)
with open(wiki_file, 'r', encoding='utf-8') as f:
    first_line = f.readline().rstrip('\n')
parts = first_line.split('\t')
print('Number of fields in first row:', len(parts))
print('Pathway field 1:', parts[0])
print('Pathway field 2:', parts[1] if len(parts) > 1 else None)
print('First 10 genes:', parts[2:12])
import json
import pandas as pd
from pathlib import Path
RAW_DIR = Path.cwd() / 'data_raw'
PROCESSED_DIR = Path.cwd() / 'data_processed'
PROCESSED_DIR.mkdir(exist_ok=True)
hallmark_rows = []
hallmark_dir = RAW_DIR / 'Hallmarks_Pathways'
for file in sorted(hallmark_dir.glob('*.json')):
    with open(file, 'r', encoding='utf-8') as f:
        obj = json.load(f)
    for pathway_name, pdata in obj.items():
        for gene in pdata['geneSymbols']:
            hallmark_rows.append({'Database': 'MSigDB_Hallmark', 'Pathway_ID': pathway_name, 'Pathway_Name': pathway_name, 'Gene_Symbol': gene})
hallmark_df = pd.DataFrame(hallmark_rows)
kegg_file = RAW_DIR / 'KEGG_Homo_sapiens_final_database' / 'KEGG_Homo_sapiens_pathway_gene_database.csv'
kegg_raw = pd.read_csv(kegg_file)
kegg_df = kegg_raw[['Pathway_ID', 'Pathway_Name', 'Gene_Symbol']].copy()
kegg_df.insert(0, 'Database', 'KEGG')
lc_file = RAW_DIR / 'LCPathways' / '41568_2020_240_MOESM4_ESM.csv'
lc_raw = pd.read_csv(lc_file)
lc_df = lc_raw[['name', 'genes']].copy()
lc_df['Gene_Symbol'] = lc_df['genes'].str.split('|')
lc_df = lc_df.explode('Gene_Symbol')
lc_df = lc_df.rename(columns={'name': 'Pathway_Name'})
pathway_ids = {name: f'LCP_{i + 1:03d}' for i, name in enumerate(lc_raw['name'])}
lc_df['Pathway_ID'] = lc_df['Pathway_Name'].map(pathway_ids)
lc_df['Database'] = 'LCPathways'
lc_df = lc_df[['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']]
reactome_file = RAW_DIR / 'ReactomePathways.gmt' / 'ReactomePathways.gmt'
reactome_rows = []
with open(reactome_file, 'r', encoding='utf-8') as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        pathway_name = parts[0]
        pathway_id = parts[1]
        genes = parts[2:]
        for gene in genes:
            reactome_rows.append({'Database': 'Reactome', 'Pathway_ID': pathway_id, 'Pathway_Name': pathway_name, 'Gene_Symbol': gene})
reactome_df = pd.DataFrame(reactome_rows)
dfs = {'Hallmark': hallmark_df, 'KEGG': kegg_df, 'LCPathways': lc_df, 'Reactome': reactome_df}
for name, df in dfs.items():
    df['Gene_Symbol'] = df['Gene_Symbol'].astype(str).str.strip()
    df.dropna(subset=['Pathway_Name', 'Gene_Symbol'], inplace=True)
    df.drop_duplicates(subset=['Pathway_ID', 'Gene_Symbol'], inplace=True)
    df.reset_index(drop=True, inplace=True)
for name, df in dfs.items():
    print('\n', '=' * 60)
    print(name)
    print('=' * 60)
    print('Rows:', len(df))
    print('Pathways:', df['Pathway_ID'].nunique())
    print('Unique genes:', df['Gene_Symbol'].nunique())
    print('\nExample:')
    print(df.head())
from pathlib import Path
RAW_DIR = Path.cwd() / 'data_raw'
wiki_file = RAW_DIR / 'wikipathways' / 'wikipathways-20230210-gmt-Homo_sapiens.gmt'
wiki_rows = []
wiki_ids = set()
wiki_pathway_count = 0
with open(wiki_file, 'r', encoding='utf-8') as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        pathway_name = parts[0]
        pathway_url = parts[1]
        gene_ids = parts[2:]
        wiki_pathway_count += 1
        for gene_id in gene_ids:
            gene_id = gene_id.strip()
            if gene_id:
                wiki_ids.add(gene_id)
print('WikiPathways pathways:', wiki_pathway_count)
print('Unique numeric gene IDs:', len(wiki_ids))
print('\nFirst 20 IDs:')
print(sorted(list(wiki_ids))[:20])
url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
import pandas as pd
from pathlib import Path
import urllib.request
RAW_DIR = Path.cwd() / 'data_raw'
gene_info_gz = RAW_DIR / 'Homo_sapiens.gene_info.gz'
url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
if not gene_info_gz.exists():
    print('Downloading NCBI Homo sapiens gene_info...')
    urllib.request.urlretrieve(url, gene_info_gz)
    print('Download complete.')
else:
    print('NCBI gene_info already exists.')
print('Saved at:', gene_info_gz)
import pandas as pd
gene_info = pd.read_csv(gene_info_gz, sep='\t', compression='gzip', dtype=str)
print('NCBI gene_info shape:', gene_info.shape)
print('Columns:', list(gene_info.columns))
ncbi_map = gene_info[['GeneID', 'Symbol']].dropna().drop_duplicates(subset=['GeneID'])
print('Unique NCBI Gene IDs:', ncbi_map['GeneID'].nunique())
print('\nExample:')
print(ncbi_map.head())
from pathlib import Path
import pandas as pd
RAW_DIR = Path.cwd() / 'data_raw'
wiki_file = RAW_DIR / 'wikipathways' / 'wikipathways-20230210-gmt-Homo_sapiens.gmt'
wiki_rows = []
with open(wiki_file, 'r', encoding='utf-8') as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        field1 = parts[0]
        gene_ids = parts[2:]
        info = field1.split('%')
        pathway_name = info[0]
        pathway_id = None
        for x in info:
            if x.startswith('WP'):
                pathway_id = x
                break
        for gene_id in gene_ids:
            gene_id = gene_id.strip()
            if gene_id:
                wiki_rows.append({'Database': 'WikiPathways', 'Pathway_ID': pathway_id, 'Pathway_Name': pathway_name, 'NCBI_Gene_ID': gene_id})
wiki_raw_df = pd.DataFrame(wiki_rows)
print('Raw WikiPathways rows:', len(wiki_raw_df))
print('Raw WikiPathways pathways:', wiki_raw_df['Pathway_ID'].nunique())
print('Unique WikiPathways gene IDs:', wiki_raw_df['NCBI_Gene_ID'].nunique())
wiki_mapped = wiki_raw_df.merge(ncbi_map, left_on='NCBI_Gene_ID', right_on='GeneID', how='left')
wiki_mapped = wiki_mapped.rename(columns={'Symbol': 'Gene_Symbol'})
total_unique_ids = wiki_mapped['NCBI_Gene_ID'].nunique()
mapped_unique_ids = wiki_mapped.loc[wiki_mapped['Gene_Symbol'].notna(), 'NCBI_Gene_ID'].nunique()
unmapped_unique_ids = total_unique_ids - mapped_unique_ids
mapping_rate = mapped_unique_ids / total_unique_ids * 100
print('\n' + '=' * 60)
print('WIKIPATHWAYS MAPPING REPORT')
print('=' * 60)
print('Total unique IDs:', total_unique_ids)
print('Mapped unique IDs:', mapped_unique_ids)
print('Unmapped unique IDs:', unmapped_unique_ids)
print(f'Mapping rate: {mapping_rate:.2f}%')
wiki_df = wiki_mapped[['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']].copy()
wiki_df = wiki_df.dropna(subset=['Pathway_ID', 'Pathway_Name', 'Gene_Symbol'])
wiki_df['Gene_Symbol'] = wiki_df['Gene_Symbol'].astype(str).str.strip()
wiki_df = wiki_df.drop_duplicates(subset=['Pathway_ID', 'Gene_Symbol'])
wiki_df = wiki_df.reset_index(drop=True)
print('\nFinal WikiPathways rows:', len(wiki_df))
print('Final WikiPathways pathways:', wiki_df['Pathway_ID'].nunique())
print('Final WikiPathways unique genes:', wiki_df['Gene_Symbol'].nunique())
print('\nExample:')
print(wiki_df.head())
PROCESSED_DIR = Path.cwd() / 'data_processed'
PROCESSED_DIR.mkdir(exist_ok=True)
hallmark_df.to_csv(PROCESSED_DIR / 'Hallmark_standardized.csv', index=False)
kegg_df.to_csv(PROCESSED_DIR / 'KEGG_standardized.csv', index=False)
lc_df.to_csv(PROCESSED_DIR / 'LCPathways_standardized.csv', index=False)
reactome_df.to_csv(PROCESSED_DIR / 'Reactome_standardized.csv', index=False)
wiki_df.to_csv(PROCESSED_DIR / 'WikiPathways_standardized.csv', index=False)
print('Saved 5 standardized files.')
master_df = pd.concat([hallmark_df, kegg_df, lc_df, reactome_df, wiki_df], ignore_index=True)
master_df = master_df.drop_duplicates(subset=['Database', 'Pathway_ID', 'Gene_Symbol']).reset_index(drop=True)
master_file = PROCESSED_DIR / 'all_pathway_resources_master.csv'
master_df.to_csv(master_file, index=False)
print('Master file saved:', master_file)
print('Master shape:', master_df.shape)
print('\nPathways by database:')
print(master_df.groupby('Database')['Pathway_ID'].nunique().sort_values())
print('\nUnique genes by database:')
print(master_df.groupby('Database')['Gene_Symbol'].nunique().sort_values())
pathway_sizes = master_df.groupby(['Database', 'Pathway_ID', 'Pathway_Name'])['Gene_Symbol'].nunique().reset_index(name='Gene_Count')
print('Total pathways:', len(pathway_sizes))
print('\nExample:')
print(pathway_sizes.head())
summary_stats = pathway_sizes.groupby('Database').agg(Number_of_Pathways=('Pathway_ID', 'nunique'), Mean_Genes_per_Pathway=('Gene_Count', 'mean'), Median_Genes_per_Pathway=('Gene_Count', 'median'), Min_Genes_per_Pathway=('Gene_Count', 'min'), Max_Genes_per_Pathway=('Gene_Count', 'max')).reset_index()
unique_genes = master_df.groupby('Database')['Gene_Symbol'].nunique().reset_index(name='Unique_Genes')
summary_stats = summary_stats.merge(unique_genes, on='Database', how='left')
summary_stats['Mean_Genes_per_Pathway'] = summary_stats['Mean_Genes_per_Pathway'].round(2)
print(summary_stats)
summary_stats.to_csv(PROCESSED_DIR / 'database_structural_summary.csv', index=False)
pathway_sizes.to_csv(PROCESSED_DIR / 'pathway_size_by_database.csv', index=False)
print('Saved structural summary files.')
from pathlib import Path
import pandas as pd
from datetime import date
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / 'data_processed'
RESULTS_DIR = PROJECT_ROOT / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
print('Project root:')
print(PROJECT_ROOT)
print()
resource_files = {'LCPathways': PROCESSED_DIR / 'LCPathways_standardized.csv', 'KEGG': PROCESSED_DIR / 'KEGG_standardized.csv', 'Reactome': PROCESSED_DIR / 'Reactome_standardized.csv', 'WikiPathways': PROCESSED_DIR / 'WikiPathways_standardized.csv', 'MSigDB Hallmark': PROCESSED_DIR / 'Hallmark_standardized.csv'}
print('Checking standardized files...\n')
missing_files = []
for resource, path in resource_files.items():
    if path.exists():
        print(f'✓ {resource:<18} {path.name}')
    else:
        print(f'✗ {resource:<18} FILE NOT FOUND: {path}')
        missing_files.append(str(path))
if missing_files:
    raise FileNotFoundError('\nOne or more standardized files are missing.\nDo not continue until the paths are corrected.')
print('\nAll 5 standardized files were found.')
required_columns = ['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']
resource_data = {}
summary_rows = []
print('\n' + '=' * 75)
print('VALIDATING STANDARDIZED RESOURCE FILES')
print('=' * 75)
for resource, path in resource_files.items():
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f'{resource} is missing required columns: {missing_columns}')
    df = df[required_columns].copy()
    for col in required_columns:
        df[col] = df[col].astype(str).str.strip()
    df['Gene_Symbol'] = df['Gene_Symbol'].str.upper()
    invalid_values = {'', 'NAN', 'NONE', 'NA'}
    for col in required_columns:
        df = df[~df[col].str.upper().isin(invalid_values)]
    before = len(df)
    df = df.drop_duplicates(subset=['Database', 'Pathway_ID', 'Gene_Symbol']).reset_index(drop=True)
    duplicates_removed = before - len(df)
    resource_data[resource] = df
    n_pathways = df['Pathway_ID'].nunique()
    n_genes = df['Gene_Symbol'].nunique()
    n_rows = len(df)
    summary_rows.append({'Resource': resource, 'Pathway_Count': n_pathways, 'Unique_Gene_Count': n_genes, 'Pathway_Gene_Associations': n_rows, 'Duplicates_Removed_During_Validation': duplicates_removed})
    print(f'\n{resource}')
    print('-' * len(resource))
    print(f'Rows:          {n_rows:,}')
    print(f'Pathways:      {n_pathways:,}')
    print(f'Unique genes:  {n_genes:,}')
    print(f'Duplicates removed: {duplicates_removed:,}')
resource_metadata = pd.DataFrame([{'Database': 'LCPathways', 'Full_Name': 'LCPathways', 'Version_or_Release': 'Supplementary dataset from 2020 publication', 'Current_Study_File': '41568_2020_240_MOESM4_ESM.csv', 'Source': 'LCPathways publication supplementary resource', 'Gene_Identifier_Used': 'Gene Symbol', 'Resource_Type': 'Cancer-focused pathway resource', 'Scope': 'Cancer-focused', 'Hierarchical': 'No explicit hierarchy used in current analysis', 'Notes': '241 pathways in the downloaded dataset. Original pathway names and source information retained.'}, {'Database': 'KEGG', 'Full_Name': 'Kyoto Encyclopedia of Genes and Genomes (KEGG)', 'Version_or_Release': 'Retrieved 2026-08-15', 'Current_Study_File': 'KEGG Homo sapiens REST-derived dataset', 'Source': 'KEGG REST API', 'Gene_Identifier_Used': 'Gene Symbol', 'Resource_Type': 'Curated pathway database', 'Scope': 'Broad', 'Hierarchical': 'Pathways organized into KEGG functional categories', 'Notes': 'Human pathways retrieved from KEGG and converted to standardized pathway-gene associations.'}, {'Database': 'Reactome', 'Full_Name': 'Reactome', 'Version_or_Release': 'Release not yet recorded from downloaded GMT', 'Current_Study_File': 'ReactomePathways.gmt', 'Source': 'Reactome pathway GMT distribution', 'Gene_Identifier_Used': 'Gene Symbol', 'Resource_Type': 'Curated pathway database', 'Scope': 'Broad', 'Hierarchical': 'Yes', 'Notes': 'Reactome has a strong hierarchical pathway organization. Current GMT contains pathways at multiple levels of granularity. Exact Reactome release should be recorded before final manuscript.'}, {'Database': 'WikiPathways', 'Full_Name': 'WikiPathways', 'Version_or_Release': '2023-02-10', 'Current_Study_File': 'wikipathways-20230210-gmt-Homo_sapiens.gmt', 'Source': 'WikiPathways Homo sapiens GMT distribution', 'Gene_Identifier_Used': 'NCBI Gene ID in raw GMT; mapped to Gene Symbol', 'Resource_Type': 'Community-curated pathway database', 'Scope': 'Broad', 'Hierarchical': 'No single comprehensive hierarchy used', 'Notes': '7,985 unique mapped genes. 7,985 of 7,990 raw NCBI Gene IDs were mapped (99.94%). This 2023 release should be updated before the final BRCA-PRISM analysis.'}, {'Database': 'MSigDB Hallmark', 'Full_Name': 'MSigDB Hallmark Gene Sets', 'Version_or_Release': '2025.1.Hs', 'Current_Study_File': 'Hallmarks_Pathways/*.json', 'Source': 'Molecular Signatures Database (MSigDB)', 'Gene_Identifier_Used': 'Gene Symbol', 'Resource_Type': 'Hallmark gene-set collection', 'Scope': 'Broad biological processes / cancer-relevant hallmarks', 'Hierarchical': 'No', 'Notes': 'Hallmark is a gene-set collection rather than a traditional pathway database. Current downloaded files produced 49 unique Hallmark sets because HALLMARK_OXIDATIVE_PHOSPHORYLATION appeared twice. Resolve to the official 50 sets before the final analysis.'}])
resource_summary = pd.DataFrame(summary_rows)
resource_summary = resource_summary.sort_values('Resource').reset_index(drop=True)
step1_resource_table = resource_metadata.merge(resource_summary, left_on='Database', right_on='Resource', how='left')
step1_resource_table = step1_resource_table.drop(columns=['Resource'])
master = pd.concat(list(resource_data.values()), ignore_index=True)
master = master.drop_duplicates(subset=['Database', 'Pathway_ID', 'Gene_Symbol']).reset_index(drop=True)
master = master[['Database', 'Pathway_ID', 'Pathway_Name', 'Gene_Symbol']]
for resource, df in resource_data.items():
    safe_name = resource.replace(' ', '_').replace('/', '_')
    output_file = PROCESSED_DIR / f'{safe_name}_standardized_validated.csv'
    df.to_csv(output_file, index=False)
master_file = PROCESSED_DIR / 'BRCA_PRISM_all_pathway_resources_master.csv'
metadata_file = RESULTS_DIR / 'BRCA_PRISM_resource_metadata.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_resource_summary.csv'
step1_table_file = RESULTS_DIR / 'BRCA_PRISM_step1_resource_documentation.csv'
master.to_csv(master_file, index=False)
resource_metadata.to_csv(metadata_file, index=False)
resource_summary.to_csv(summary_file, index=False)
step1_resource_table.to_csv(step1_table_file, index=False)
print('\n' + '=' * 75)
print('STEP 1 FINAL RESOURCE SUMMARY')
print('=' * 75)
display(step1_resource_table[['Database', 'Version_or_Release', 'Pathway_Count', 'Unique_Gene_Count', 'Pathway_Gene_Associations', 'Gene_Identifier_Used', 'Scope', 'Hierarchical']])
print('\nMaster standardized table:')
print('Shape:', master.shape)
print('\nColumns:')
print(master.columns.tolist())
print('\nRows by resource:')
print(master['Database'].value_counts())
print('\nNumber of pathways by resource:')
print(master.groupby('Database')['Pathway_ID'].nunique().sort_values(ascending=False))
print('\nNumber of unique genes by resource:')
print(master.groupby('Database')['Gene_Symbol'].nunique().sort_values(ascending=False))
print('\nSaved:')
print(master_file)
print(metadata_file)
print(summary_file)
print(step1_table_file)
print('\nSTEP 1 COMPLETE.')
from pathlib import Path
import pandas as pd
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / 'data_processed'
RESULTS_DIR = PROJECT_ROOT / 'results'
master_file = PROCESSED_DIR / 'BRCA_PRISM_all_pathway_resources_master.csv'
metadata_file = RESULTS_DIR / 'BRCA_PRISM_resource_metadata.csv'
summary_file = RESULTS_DIR / 'BRCA_PRISM_resource_summary.csv'
documentation_file = RESULTS_DIR / 'BRCA_PRISM_step1_resource_documentation.csv'
master = pd.read_csv(master_file)
metadata = pd.read_csv(metadata_file)
summary = pd.read_csv(summary_file)
documentation = pd.read_csv(documentation_file)
name_map = {'MSigDB_Hallmark': 'MSigDB Hallmark', 'Hallmark': 'MSigDB Hallmark'}
master['Database'] = master['Database'].replace(name_map)
if 'Database' in metadata.columns:
    metadata['Database'] = metadata['Database'].replace(name_map)
if 'Resource' in summary.columns:
    summary['Resource'] = summary['Resource'].replace(name_map)
if 'Database' in documentation.columns:
    documentation['Database'] = documentation['Database'].replace(name_map)
master.to_csv(master_file, index=False)
metadata.to_csv(metadata_file, index=False)
summary.to_csv(summary_file, index=False)
documentation.to_csv(documentation_file, index=False)
print('Databases in master:')
print(master['Database'].unique())
print('\nFinal counts:')
final_check = master.groupby('Database').agg(Pathways=('Pathway_ID', 'nunique'), Unique_Genes=('Gene_Symbol', 'nunique'), Associations=('Gene_Symbol', 'size')).reset_index()
display(final_check)
print('\nMaster shape:', master.shape)
print('\nSTEP 1 is now frozen for the current working dataset.')
