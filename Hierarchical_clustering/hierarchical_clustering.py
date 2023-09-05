import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import matplotlib as mpl
import pandas as pd
import seaborn as sns
##import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

algo_list=['KMC','Gusev','Annable','Eichinger', 'Grimson_new', 'Lei', 'DWS_1', 'DWS_2', 'DWS_3', 'DWS_4']#, 'DWS_2', 'DWS_3', 'DWS_4']
num_algo=len(algo_list)
num_runs=3
EMD_mean_mat=np.zeros((num_algo,num_algo))
EMD_std_mat=np.zeros((num_algo,num_algo))
for n1 in range(0,num_algo):
    for n2 in range(n1+1,num_algo):
        net1=algo_list[n1]
        net2=algo_list[n2]
##        num_replicates=0
        EMD_replicates=[]
        for run1 in range(1,num_runs+1):
            for run2 in range(1,num_runs+1):
                filename='EMD_mat_'+net1+'_'+str(run1)+'_'+net2+'_'+str(run2)+'.txt'
##                num_replicates=num_replicates+1
                try:
                    emd=np.genfromtxt(filename)
                except FileNotFoundError:
                    filename='EMD_mat_'+net2+'_'+str(run2)+'_'+net1+'_'+str(run1)+'.txt'
                    emd=np.genfromtxt(filename)
                emd=emd[0,1]
                EMD_replicates.append(emd)
        EMD_mean=np.mean(np.array(EMD_replicates))
        EMD_std=np.std(np.array(EMD_replicates))
        EMD_mean_mat[n1,n2]=EMD_mean
        EMD_mean_mat[n2,n1]=EMD_mean
        EMD_std_mat[n1,n2]=EMD_std
        EMD_std_mat[n2,n1]=EMD_std

np.savetxt('EMD_mean_mat_with_DWS.txt',EMD_mean_mat)
np.savetxt('EMD_std_mat_with_DWS.txt',EMD_std_mat)

algo=['KMC','Gusev','Annable','Leung', 'Grimson', 'Lei', 'DWS_1', 'DWS_2', 'DWS_3', 'DWS_4']#, 'DWS_2', 'DWS_3', 'DWS_4']




similarity=EMD_mean_mat




D=similarity

condensedD=[]#np.zeros(comb(len(D),2))
for i in range(0,len(D)):
    for j in range(i+1,len(D)):
        condensedD.append(D[i,j])


# Compute and plot first dendrogram.

method='ward'
fig, ax = plt.subplots()
Y = sch.linkage(condensedD, method=method, optimal_ordering=True)

for label in (ax.get_xticklabels() ):
    label.set_fontsize(20)
    label.set_fontweight('bold')
##    label.set_rotation(45)
ax.tick_params(width=2)
for label in (ax.get_yticklabels() ):
    label.set_fontsize(15)
    label.set_fontweight('bold')

[x.set_linewidth(1.5) for x in ax.spines.values()]

with plt.rc_context({'lines.linewidth': 3.5}):
    Z1 = sch.dendrogram(Y,p=30,color_threshold= 0.7*max(D[:,2]),labels=algo, ax=ax)#, orientation='left')


ax.set_ylabel('Distance between clusters', fontsize=15, fontweight='bold')
plt.show()


plt.figure()
# Compute and plot second dendrogram.

Y = sch.linkage(condensedD, method=method, optimal_ordering=True)
Z2 = sch.dendrogram(Y, orientation='right')


# Plot distance matrix.

fig,axmatrix=plt.subplots()

idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1,:]
D = D[:,idx2]
algo=[algo[i] for i in idx1]
print(algo)
similarity = pd.DataFrame(D, index = algo, columns = algo)
im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
axmatrix.set_xticks([])  # remove axis labels
axmatrix.set_yticks([])  # remove axis labels
axmatrix.set_xticks(axmatrix.get_xticks())
axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=45, ha='right')
axmatrix.set_yticks(axmatrix.get_yticks())
axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=45)


fig, ax = plt.subplots()#figsize=(13,10))

# Add title to the Heat map
title = "Similarity of algorithms"

# Set the font size and the distance of the title from the plot
plt.title(title,fontsize=18)

# Use the heatmap function from the seaborn package
sns.heatmap(similarity,annot=True,cmap='GnBu',linewidths=0.30, linecolor='black',annot_kws={"fontsize":12, "weight":20},ax=ax)#,xticklabels=x_axis_labels, yticklabels=y_axis_labels)
plt.savefig('similarity_algo')

# Display the Pharma Sector Heatmap
plt.show()

