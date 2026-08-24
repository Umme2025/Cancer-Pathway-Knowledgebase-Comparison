from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / 'results'
coverage = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step8_TCGA_family_coverage.csv')
pam50 = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step9B_PAM50_family_association.csv')
stability = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step10B_gene_overlap_vs_patient_stability.csv')
patient_exp = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step14B_patient_explanation_pairwise_agreement.csv')
ablation = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step16C_repeatedCV_ablation_summary.csv')
ablation_stats = pd.read_csv(RESULTS_DIR / 'BRCA_PRISM_step16D_paired_ablation_statistics.csv')
rho, p = spearmanr(stability['Gene_Jaccard'], stability['Spearman_Rho'])
identical_top3 = (patient_exp['Shared_Top3_Families'] == 3).mean()
top_family = pam50.sort_values('Epsilon_Squared', ascending=False).iloc[0]
macro = dict(zip(ablation['Representation'], ablation['Mean_Macro_F1']))
rows = [{'Finding': 'TCGA coverage of family union genes', 'Metric': 'Mean retention', 'Value': coverage['Union_Retention'].mean(), 'Interpretation': 'BRCA-PRISM family genes are well represented in TCGA-BRCA.'}, {'Finding': 'TCGA coverage of consensus genes', 'Metric': 'Mean retention', 'Value': coverage['Consensus_Retention'].mean(), 'Interpretation': 'Cross-resource consensus genes are almost completely measurable.'}, {'Finding': 'Strongest PAM50-associated family', 'Metric': f"{top_family['Family_ID']} epsilon-squared", 'Value': top_family['Epsilon_Squared'], 'Interpretation': 'DNA-replication-related activity showed the strongest subtype association.'}, {'Finding': 'Gene-set similarity vs patient-level stability', 'Metric': 'Spearman rho', 'Value': rho, 'Interpretation': 'More similar pathway gene definitions produce more stable patient-level activity.'}, {'Finding': 'Gene-set similarity vs patient-level stability', 'Metric': 'P-value', 'Value': p, 'Interpretation': 'The association is highly statistically significant.'}, {'Finding': 'Patient-level explanation stability', 'Metric': 'Fraction with identical Top-3 families', 'Value': identical_top3, 'Interpretation': 'Most patient explanations change at least partially when the pathway resource changes.'}, {'Finding': 'Repeated-CV ablation', 'Metric': 'Resource-specific genes Macro-F1', 'Value': macro['Resource_Specific_Genes'], 'Interpretation': 'Resource-specific pathway content retained the strongest subtype-discriminative signal.'}, {'Finding': 'Repeated-CV ablation', 'Metric': 'All family genes Macro-F1', 'Value': macro['All_Family_Genes'], 'Interpretation': 'Using all family genes produced intermediate predictive performance.'}, {'Finding': 'Repeated-CV ablation', 'Metric': 'Consensus genes Macro-F1', 'Value': macro['Consensus_Genes'], 'Interpretation': 'Consensus genes provided stable shared biology but lower subtype discrimination.'}]
summary = pd.DataFrame(rows)
for _, r in ablation_stats.iterrows():
    summary.loc[len(summary)] = {'Finding': 'Paired ablation comparison', 'Metric': r['Comparison'], 'Value': r['FDR'], 'Interpretation': 'FDR-adjusted Wilcoxon test across 50 paired CV splits.'}
print('=' * 100)
print('BRCA-PRISM FINAL MAIN RESULTS SUMMARY')
print('=' * 100)
display(summary)
OUT_FILE = RESULTS_DIR / 'BRCA_PRISM_FINAL_main_results_summary.csv'
summary.to_csv(OUT_FILE, index=False)
print('\nMAIN RESULTS FROZEN')
print('Saved:', OUT_FILE)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
FIGURE_DIR = Path('C:\\Users\\aps211\\OneDrive - The University of Texas-Rio Grande Valley\\Vidio Presentation\\Research_Own\\Pathway Comparison\\figures')
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
all_patients = pd.DataFrame({'Pair': ['KEGG vs\nConsensus', 'KEGG vs\nReactome', 'KEGG vs\nWikiPathways', 'Reactome vs\nConsensus', 'Reactome vs\nWikiPathways', 'WikiPathways vs\nConsensus'], 'Prediction_Agreement': [0.7813, 0.6078, 0.7146, 0.6376, 0.6253, 0.807], 'Mean_Shared_Top3': [2.1869, 1.7977, 1.8747, 1.8809, 1.7248, 2.0277], 'Mean_Top3_Jaccard': [0.6264, 0.4774, 0.4908, 0.5075, 0.443, 0.5486], 'Mean_Contribution_Rank_Rho': [0.6792, 0.4171, 0.502, 0.5145, 0.4071, 0.604]})
same_prediction = pd.DataFrame({'Pair': ['KEGG vs\nConsensus', 'KEGG vs\nReactome', 'KEGG vs\nWikiPathways', 'Reactome vs\nConsensus', 'Reactome vs\nWikiPathways', 'WikiPathways vs\nConsensus'], 'N_Same_Prediction': [761, 592, 696, 621, 609, 786], 'Mean_Shared_Top3': [2.3443, 2.0439, 1.9353, 2.1224, 1.8818, 2.1005], 'Mean_Top3_Jaccard': [0.6853, 0.5625, 0.5111, 0.5902, 0.4943, 0.5766], 'Mean_Contribution_Rank_Rho': [0.7511, 0.5553, 0.5272, 0.6612, 0.5056, 0.646]})
all_patients.to_csv(FIGURE_DIR / 'BRCA_PRISM_step14B_all_patients_summary.csv', index=False)
same_prediction.to_csv(FIGURE_DIR / 'BRCA_PRISM_step14B_same_prediction_summary.csv', index=False)
fig, axes = plt.subplots(2, 1, figsize=(14, 10), constrained_layout=True)
ax = axes[0]
x = np.arange(len(all_patients))
width = 0.24
metrics_A = [('Prediction_Agreement', 'Prediction agreement'), ('Mean_Top3_Jaccard', 'Top-3 Jaccard'), ('Mean_Contribution_Rank_Rho', 'Contribution-rank rho')]
for i, (col, label) in enumerate(metrics_A):
    ax.bar(x + (i - 1) * width, all_patients[col], width=width, label=label)
ax.set_title('A. Explanation stability across all TCGA-BRCA patients', fontsize=14, weight='bold')
ax.set_ylabel('Value')
ax.set_xticks(x)
ax.set_xticklabels(all_patients['Pair'])
ax.set_ylim(0, 1.0)
ax.legend(frameon=False, ncol=3)
ax.grid(axis='y', linestyle='--', alpha=0.4)
ax = axes[1]
x2 = np.arange(len(same_prediction))
width2 = 0.28
metrics_B = [('Mean_Top3_Jaccard', 'Top-3 Jaccard'), ('Mean_Contribution_Rank_Rho', 'Contribution-rank rho')]
for i, (col, label) in enumerate(metrics_B):
    ax.bar(x2 + (i - 0.5) * width2, same_prediction[col], width=width2, label=label)
ax.set_title('B. Explanation stability among patients with the same predicted PAM50 subtype', fontsize=14, weight='bold')
ax.set_ylabel('Value')
ax.set_xticks(x2)
ax.set_xticklabels(same_prediction['Pair'])
ax.set_ylim(0, 1.0)
ax.legend(frameon=False, ncol=2)
ax.grid(axis='y', linestyle='--', alpha=0.4)
out_png = FIGURE_DIR / 'BRCA_PRISM_Figure3_patient_explanation_stability.png'
out_pdf = FIGURE_DIR / 'BRCA_PRISM_Figure3_patient_explanation_stability.pdf'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.show()
print('Saved figure:')
print(out_png)
print(out_pdf)
print('\nSummary for patients with SAME predicted subtype:')
display(same_prediction[['Pair', 'N_Same_Prediction', 'Mean_Shared_Top3', 'Mean_Top3_Jaccard', 'Mean_Contribution_Rank_Rho']])
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
perf_df = pd.DataFrame({'Representation': ['Resource-Specific\nGenes', 'All Family\nGenes', 'Consensus\nGenes'], 'Accuracy': [0.6908, 0.6649, 0.6596], 'Balanced Accuracy': [0.6788, 0.636, 0.6237], 'Macro-F1': [0.5988, 0.5634, 0.5553], 'SD_Accuracy': [0.034, 0.0323, 0.0288], 'SD_Balanced_Accuracy': [0.0493, 0.0536, 0.0461], 'SD_Macro_F1': [0.0384, 0.0379, 0.0313]})
paired_df = pd.DataFrame({'Comparison': ['Resource-Specific vs All Family', 'Resource-Specific vs Consensus', 'All Family vs Consensus'], 'Mean_Paired_Difference': [0.035369, 0.043514, 0.008145], 'FDR': [2.427494e-08, 2.427494e-08, 0.08054833]})
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={'width_ratios': [1.05, 1.1]})
ax = axes[0]
x = np.arange(len(perf_df))
width = 0.22
bars1 = ax.bar(x - width, perf_df['Accuracy'], width, yerr=perf_df['SD_Accuracy'], capsize=4, label='Accuracy')
bars2 = ax.bar(x, perf_df['Balanced Accuracy'], width, yerr=perf_df['SD_Balanced_Accuracy'], capsize=4, label='Balanced Accuracy')
bars3 = ax.bar(x + width, perf_df['Macro-F1'], width, yerr=perf_df['SD_Macro_F1'], capsize=4, label='Macro-F1')
ax.set_title('A. Repeated-CV Classification Performance', fontsize=13, fontweight='bold')
ax.set_ylabel('Performance', fontsize=11)
ax.set_ylim(0, 0.82)
ax.set_xticks(x)
ax.set_xticklabels(perf_df['Representation'], fontsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.grid(axis='y', linestyle='--', alpha=0.35)
ax.set_axisbelow(True)
bar_sets = [(bars1, perf_df['SD_Accuracy']), (bars2, perf_df['SD_Balanced_Accuracy']), (bars3, perf_df['SD_Macro_F1'])]
for bars, errors in bar_sets:
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + errors.iloc[i] + 0.012, f'{height:.3f}', ha='center', va='bottom', fontsize=9)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=10)
ax = axes[1]
y = np.arange(len(paired_df))
bars = ax.barh(y, paired_df['Mean_Paired_Difference'], height=0.55)
ax.set_title('B. Paired Macro-F1 Comparison', fontsize=13, fontweight='bold')
ax.set_xlabel('Mean Paired Macro-F1 Difference', fontsize=11)
ax.set_yticks(y)
ax.set_yticklabels(paired_df['Comparison'], fontsize=10)
ax.tick_params(axis='x', labelsize=10)
ax.grid(axis='x', linestyle='--', alpha=0.35)
ax.set_axisbelow(True)
ax.set_xlim(0, 0.07)
ax.invert_yaxis()
for i, row in paired_df.iterrows():
    diff = row['Mean_Paired_Difference']
    fdr = row['FDR']
    if fdr < 0.001:
        fdr_text = 'FDR = $2.43\\times10^{-8}$'
    else:
        fdr_text = f'FDR = {fdr:.4f}'
    significance_text = 'Significant' if fdr < 0.05 else 'Not significant'
    ax.text(diff + 0.0013, i, f'Δ = {diff:.4f}\n{fdr_text}\n{significance_text}', va='center', ha='left', fontsize=8.7, linespacing=1.12)
plt.tight_layout(rect=[0, 0.06, 1, 1], w_pad=3.5)
PNG_FILE = FIGURES_DIR / 'BRCA_PRISM_Figure4_gene_ablation_improved.png'
PDF_FILE = FIGURES_DIR / 'BRCA_PRISM_Figure4_gene_ablation_improved.pdf'
plt.savefig(PNG_FILE, dpi=300, bbox_inches='tight', pad_inches=0.15)
plt.savefig(PDF_FILE, bbox_inches='tight', pad_inches=0.15)
plt.show()
print('Saved:')
print(PNG_FILE)
print(PDF_FILE)
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = PROJECT_ROOT / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)

def rounded_box(ax, x, y, w, h, text, facecolor, fontsize=11, fontweight='normal', edgecolor='0.25', linewidth=1.4):
    box = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.015,rounding_size=0.025', facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth)
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fontsize, fontweight=fontweight, linespacing=1.25)
    return box

def arrow(ax, x1, y1, x2, y2):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=17, linewidth=1.6, color='0.25', connectionstyle='arc3')
    ax.add_patch(a)
fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
ax.text(0.5, 0.965, 'BRCA-PRISM: Cross-Resource Pathway Integration and Stability Analysis', ha='center', va='center', fontsize=20, fontweight='bold')
ax.text(0.5, 0.895, '1. Pathway / Gene-Set Resources', ha='center', fontsize=13, fontweight='bold')
resource_names = ['LCPathways', 'KEGG', 'Reactome', 'WikiPathways', 'MSigDB\nHallmark']
resource_x = [0.055, 0.245, 0.435, 0.625, 0.815]
resource_colors = ['#FCE8E6', '#E8F1F8', '#E5F4EF', '#FFF4D6', '#F1E8F5']
for x, name, color in zip(resource_x, resource_names, resource_colors):
    rounded_box(ax=ax, x=x, y=0.81, w=0.13, h=0.06, text=name, facecolor=color, fontsize=11, fontweight='bold')
for x in resource_x:
    arrow(ax, x + 0.065, 0.81, 0.5, 0.755)
rounded_box(ax=ax, x=0.27, y=0.675, w=0.46, h=0.075, text='2. Cross-Resource Harmonization and Mapping\nGene harmonization  •  pathway similarity  •  relationship screening', facecolor='#F5F5F5', fontsize=11, fontweight='bold')
arrow(ax, 0.5, 0.675, 0.5, 0.615)
rounded_box(ax=ax, x=0.25, y=0.525, w=0.5, h=0.085, text='3. BRCA-PRISM Pathway-Family Layer\n31 cross-resource pathway families', facecolor='#DCECF8', fontsize=12, fontweight='bold')
gene_boxes = [('Union\nGenes', 0.17), ('Consensus\nGenes', 0.35), ('Core\nGenes', 0.53), ('Resource-Specific\nGenes', 0.71)]
for text, x in gene_boxes:
    rounded_box(ax=ax, x=x, y=0.405, w=0.13, h=0.065, text=text, facecolor='#EAF6F2', fontsize=10, fontweight='bold')
    arrow(ax, 0.5, 0.525, x + 0.065, 0.47)
rounded_box(ax=ax, x=0.34, y=0.285, w=0.32, h=0.075, text='4. TCGA-BRCA Multi-Omics\nRNA + CNV  |  974 patients  |  5 PAM50 subtypes', facecolor='#FCE9D3', fontsize=11, fontweight='bold')
for _, x in gene_boxes:
    arrow(ax, x + 0.065, 0.405, 0.5, 0.36)
ax.text(0.5, 0.225, '5. Downstream Evaluation', ha='center', fontsize=13, fontweight='bold')
output_boxes = [('PAM50\nAssociation', 0.055, '#F9E0E0'), ('Prediction\nStability', 0.245, '#E9F3E7'), ('Patient-Level\nInterpretation', 0.435, '#FFF0D8'), ('Gene-Level\nInterpretation', 0.625, '#EAEAF8'), ('Consensus vs\nResource-Specific', 0.815, '#F2E8F7')]
for text, x, color in output_boxes:
    rounded_box(ax=ax, x=x, y=0.11, w=0.13, h=0.075, text=text, facecolor=color, fontsize=10, fontweight='bold')
    arrow(ax, 0.5, 0.285, x + 0.065, 0.185)
ax.text(0.5, 0.045, 'Goal: quantify how pathway-resource choice affects prediction and biological interpretation in breast cancer.', ha='center', fontsize=11, fontstyle='italic')
PNG_FILE = FIGURES_DIR / 'BRCA_PRISM_Figure1_workflow_clean.png'
PDF_FILE = FIGURES_DIR / 'BRCA_PRISM_Figure1_workflow_clean.pdf'
plt.savefig(PNG_FILE, dpi=300, bbox_inches='tight', pad_inches=0.2)
plt.savefig(PDF_FILE, bbox_inches='tight', pad_inches=0.2)
plt.show()
print('Saved:')
print(PNG_FILE)
print(PDF_FILE)
