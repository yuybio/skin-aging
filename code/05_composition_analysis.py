# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =============== read metadata  ===============
metadata = pd.read_csv('penis.metadata.csv', index_col=0)

final_df = pd.DataFrame()
results_df= pd.DataFrame()
results = []
pivot_table = {}
pivot_sig = {}


for cell_type in metadata['cell_type'].cat.categories:
    print(cell_type)
    cell_type_metadata = metadata[metadata['cell_type'] == cell_type]
    counts_per_sample = cell_type_metadata.groupby(['sample_id', 'cell_type']).size().reset_index(name='sample_observed')
    counts_per_sample = counts_per_sample.merge(cell_type_metadata[['sample_id', 'stage']].drop_duplicates(), on='sample_id')
        
        
    total_observed = cell_type_metadata.groupby('cell_type').size().reset_index(name='Total_observed')
    total_observed = total_observed.merge(cell_type_metadata[['sample_id','cell_type']].drop_duplicates(), on='cell_type')
        
        
    actual_counts = pd.merge(counts_per_sample, total_observed, on=['cell_type','sample_id'])
    actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

    conut_per_sample_theoretical_counts = metadata.groupby(['sample_id', '_scvi_batch']).size().reset_index(name='sample_expected')
    conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts[conut_per_sample_theoretical_counts['_scvi_batch'] == 0]
    conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts.merge(cell_type_metadata[['sample_id','stage']].drop_duplicates(), on=['sample_id'])
        
    total_expected = metadata.groupby('_scvi_batch').size().reset_index(name='Total_expected')
    total_expected = total_expected.merge(cell_type_metadata[['stage','_scvi_batch']].drop_duplicates(), on=['_scvi_batch'])
        
    theoretical_counts = pd.merge(conut_per_sample_theoretical_counts, total_expected, on=['stage','_scvi_batch'])
       
    theoretical_counts['Percentage_expected'] = (theoretical_counts['sample_expected'] / theoretical_counts['Total_expected']) * 100
    combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','stage'])
    combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']
    combined_df=combined_df.drop(['_scvi_batch'], axis=1)
    group_stats = combined_df.groupby(['stage', 'cell_type'])['Percentage_Ratio'].agg(mean_raw='mean', std_raw='std').reset_index()
    combined_df = combined_df.merge(group_stats, on=['stage', 'cell_type'])
    combined_df_filtered = combined_df[abs(combined_df['Percentage_Ratio'] - combined_df['mean_raw']) <= 3 * combined_df['std_raw']]
    final_group_stats = combined_df_filtered.groupby(['stage', 'cell_type'])['Percentage_Ratio'].mean().reset_index().rename(columns={'Percentage_Ratio': 'Percentage_Ratio_Ave'})
    combined_df_final = combined_df_filtered.merge(final_group_stats, on=['stage', 'cell_type'])
    

    
    combined_df_final['Log_Percentage_Ratio_Ave'] = np.log2(combined_df_final['Percentage_Ratio_Ave'].replace(0, np.nan))
    combined_df_final
    combined_df.to_csv(f'NT_{cell_type}.Percentage_Ratio_raw.csv',index=False)
    combined_df_final.to_csv(f'NT_{cell_type}.Percentage_Ratio_ave.csv',index=False)
    
    
        
    for celltype_granular in combined_df_final['cell_type'].unique():
        print(celltype_granular)
        combined_df_subset = combined_df_final[combined_df_final['cell_type'] == celltype_granular]
        for target_group in combined_df_final['stage'].unique():
            print(target_group)
            target_data = combined_df_subset[combined_df_subset['stage'] == target_group]['Percentage_Ratio']
            rest_data = combined_df_subset[combined_df_subset['stage'] != target_group]['Percentage_Ratio']
            
            if not target_data.empty and not rest_data.empty:
                stat, p = mannwhitneyu(target_data, rest_data, alternative='two-sided')
                
                results.append({
                    'cell_type': cell_type,
                    'stage': target_group,
                    'Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['stage'] == target_group]['Percentage_Ratio_Ave'].values[0],
                    'Log_Percentage_Ratio_Ave':combined_df_subset[combined_df_subset['stage'] == target_group]['Log_Percentage_Ratio_Ave'].values[0],
                    'p_value': p
                })
    results_df = pd.DataFrame(results)
    corrected_p_values = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')[1]
    results_df['corrected_p_value'] = corrected_p_values
final_df = pd.concat([final_df, results_df], ignore_index=True)

final_df.to_csv('PS_ALL_cell.Percentage_Ratio_ave.csv',index=False)

pivot_table_df = final_df.pivot("cell_type", "stage", "Log_Percentage_Ratio_Ave")
pivot_sig_df = final_df.pivot("cell_type", "stage", "p_value")

Broadcelltype_order = ['KC','KC_Channel','SGC',
'FB','Pc-vSMC','VEC','LEC',
'MEL','Schwann',
'Lymphocyte','Mac-DC','LC','Mast'][::-1] 

current_order = [item for item in Broadcelltype_order if item in pivot_table_df.index]
pivot_table_df = pivot_table_df.loc[current_order]
pivot_sig_df = pivot_sig_df.loc[current_order]
pivot_table_df = pivot_table_df[['Y', 'M','O']]
pivot_sig_df = pivot_sig_df[['Y', 'M','O']]

def significance_mark(p):
    if p <= 0.05:
        return "*"
    else:
        return ""
    
sns.heatmap(pivot_table_df, annot=np.vectorize(significance_mark)(pivot_sig_df), annot_kws={"size": 10, "color": "black"}, fmt='',cmap='RdBu_r', center=0, linewidths=0.5, linecolor='white', vmin=-2, vmax=2,cbar_kws={'label': 'Log2(Ro/e)', 'shrink': 0.8})
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.savefig('PS_all_Ratio.pdf',dpi=1200,bbox_inches='tight')
plt.show()


metadata = pd.read_csv('penis.metadata.csv', index_col=0)
grouped = metadata.groupby(['sample_id', 'cell_type','stage','age']).size().reset_index(name='count')

total = grouped.groupby('sample_id')['count'].sum().reset_index(name='total_count')

merged = pd.merge(grouped, total, on='sample_id')

merged['Proportion'] = (merged['count'] / merged['total_count'] * 100).round(2)
merged['cell_type'] = merged['cell_type'].astype('category').cat.set_categories(['KC','KC_Channel','SGC', 'MEL', 'Schwann','FB', 'VEC','LEC','Pc-vSMC',
                       'Mac-DC','LC','Mast', 'Lymphocyte'])
print(merged[['sample_id', 'cell_type', 'Proportion','stage','age']])
#

df_long_sorted = merged.sort_values(by="age")
custom_colors = {
    "Y": "#FD7446FF",         
    "M": "#709AE1FF",  
    "O": "#FED439FF"    
}
g = sns.FacetGrid(df_long_sorted, col="cell_type", col_wrap=3, sharey=False, height=3.5)

# LOESS 
def draw_unified_loess(data, color, **kwargs):
    sns.regplot(
        data=data, x="age", y="Proportion", lowess=True,
        scatter=False, color="grey", line_kws={"linewidth": 2}
    )

g.map_dataframe(draw_unified_loess)

# 
def draw_colored_points(data, **kwargs):
    sns.scatterplot(
        data=data, x="age", y="Proportion", hue="stage",
        palette=custom_colors, s=40, alpha=1, legend=False
    )

g.map_dataframe(draw_colored_points)

#
g.set_titles(col_template="{col_name}")
g.set_axis_labels("Age", "Proportion (%)")
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("Cell Type Proportion vs Age (Single LOESS Fit + Colored Points)", fontsize=14)

# 
handles, labels = g.axes[0].get_legend_handles_labels()
g.fig.legend(handles, labels, title="Age Stage", loc='upper right', bbox_to_anchor=(1.05, 0.95))

# 
g.savefig("Corrected_Single_LOESS.pdf", dpi=1200)
plt.show()
