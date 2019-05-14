from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42 # Editable text when exporting pdfs
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns

umap_coordinates = pd.read_csv("results/UMAP/processed_UMAP_coordinates.csv",index_col = 0 )
predicted_cell_type = pd.read_csv("results/classification/predicted_cell_types.csv",index_col=0)
prediction_probabilities = pd.read_csv("results/classification/prediction_probabilities.csv",index_col=0)


#%% Plot prediction probabilities
plt.figure(figsize=(5,4))
plt.scatter( umap_coordinates.UMAP1, umap_coordinates.UMAP2, s = 2, 
             c = np.max( prediction_probabilities, axis = 1 ), cmap = "inferno"  )
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_ylabel("Prediction\n Probability",fontsize = 16,rotation=0)
plt.xticks([])
plt.yticks([])
plt.xlabel('UMAP1',fontsize=16)
plt.ylabel('UMAP2',fontsize=16)
plt.savefig("figures/classification/prediction_probabilities.pdf")



#%%
cell_types = np.unique( predicted_cell_type )
color_map = plt.get_cmap("tab20")
color_key = color_key = dict(zip( cell_types, color_map.colors[0:len(cell_types)] ) )
color_key["Unassigned"] = (0.8,0.8,0.8)

#%% Without filtering by predictoin confidence
plt.figure(figsize=(4,4))
for ct in cell_types:
    idx = np.squeeze(predicted_cell_type == ct)
    plt.scatter( umap_coordinates.UMAP1.values[idx],umap_coordinates.UMAP2.values[idx],s=2, c = color_key[ct])
plt.legend(cell_types,loc=(1,0),fontsize=16,markerscale = 10.)
plt.xticks([])
plt.yticks([])
plt.xlabel('UMAP1',fontsize=16)
plt.ylabel('UMAP2',fontsize=16)
plt.savefig("figures/classification/cell_type_predictions_unfiltered.pdf")

#%% After filtering 


predicted_cell_type[ np.max( prediction_probabilities, axis = 1 ) < 0.5 ] = "Unassigned"
cell_types = np.unique( predicted_cell_type )

plt.figure(figsize=(9,9))
for ct in cell_types:
    idx = np.squeeze(predicted_cell_type == ct)
    plt.scatter( umap_coordinates.UMAP1.values[idx],umap_coordinates.UMAP2.values[idx],s=2, c = color_key[ct])
plt.legend(cell_types,loc=(1,0),fontsize=16,markerscale = 10.)
plt.xticks([])
plt.yticks([])
plt.xlabel('UMAP1',fontsize=16)
plt.ylabel('UMAP2',fontsize=16)
plt.savefig("figures/classification/cell_type_predictions.pdf")

#%% Histogram


sample_type = np.array(["_".join(x.split("_")[0:3]) for x in predicted_cell_type.index.values])
unique_samples = np.unique(sample_type)


num_ct = len(cell_types)
num_samples = len(unique_samples)

ct_percentages = np.empty(( num_samples,num_ct))

for sample in unique_samples:
    sample_idx = sample_type == sample
    s = np.where(unique_samples == sample)[0]
    #print(np.sum(sample_idx))
    for ct in cell_types:
        ct_idx = predicted_cell_type.loc[:,"predicted_cell_type"].values == ct
        c = np.where( cell_types == ct)[0]
        ct_percentages[s,c] = np.sum(np.logical_and(sample_idx,ct_idx))/np.sum(sample_idx) *100

plt.figure(figsize=(4,4))    
for ct in cell_types:
    c = np.where( cell_types == ct)[0]

    if ct == cell_types[0]:
        bot = np.zeros((num_samples,))
    else:
        bot = np.squeeze(np.sum(ct_percentages[:,0:int(c)],axis=1))
    
    
    ax = plt.bar( np.arange(num_samples), np.squeeze(ct_percentages[:,c]),bottom=bot, color = color_key[ct])


plt.legend(cell_types,loc=(1,0))
plt.ylabel("Cell type percentage" )

split_labels = [("_").join( x.split("_")[0:2]) for x in unique_samples]
plt.xticks(np.arange(num_samples), split_labels, rotation = -25)
plt.savefig("figures/classification/cell_type_distribution.pdf")

df = pd.DataFrame(ct_percentages, index = unique_samples, columns=cell_types).to_csv(
        "results/classification/cell_type_percentages.csv")