import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

adata_merged = sc.read_h5ad("4Merged_samples.h5ad")

expression_data = pd.DataFrame({'CellType': adata_merged.obs['CellType'], 'RETN_Expression': adata_merged.X[:, adata_merged.var_names == 'RETN'].flatten() })

expression_data_non_zero = expression_data[expression_data['RETN_Expression'] > 0]

non_zero_cell_types = expression_data_non_zero['CellType'].unique()

adata_non_zero = adata_merged[adata_merged.obs['CellType'].isin(non_zero_cell_types)]

# Initialise figure with a size
plt.figure(figsize=(12, 6))
# Plot the violin plots without a legend and not showing the plot (it wouldn't show anyway...)
sc.pl.violin(adata_non_zero, keys='RETN', groupby='CellType', ylabel = 'RETN Expression Level', jitter = 0.4, inner = 'quartile', legend = False, show = False)
# Modify x ticks labels so they are rotated vertically, with a smaller font, and in bold
plt.xticks(rotation=90, fontsize = 5, fontweight='bold')
# Supposedly adjusts the bottom margin, but I think it's not necessary
plt.subplots_adjust(bottom = 0.3)
# Save plot at a nice resolution of 300 dpi (can be more), and frame it tight so all the information is included
plt.savefig("RETN_expression_violins.png", dpi=300, bbox_inches='tight')
plt.close()