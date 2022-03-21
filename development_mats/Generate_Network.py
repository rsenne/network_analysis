#Start by loading the functions from the referenced NetworkFunctions.py file
run NetworkFunctions.py
import pickle as pkl

#Unpickle the ROI dictionary
with open('/Users/kaitlyndorst/Documents/GitHub/networkx/Allen_Areas_dict.pickle','rb') as f:
    ROIs = pkl.load(f)
Allen_Groups = list(ROIs.values())
#ROIs = ROIs.sort_values("Allen Group Name").reset_index(drop = True)

#Then get that data
ChR2_raw_data,ChR2_nodes = loadData('/Users/kaitlyndorst/Desktop/Network_csv/ChR2_animals/ChR2_Large_Network.csv')
Control_raw_data,Control_nodes = loadData('/Users/kaitlyndorst/Desktop/Network_csv/Control_animals/Control_Large_Network.csv')

#Get the correlation and adjusted p values for data_ChR2 and data_Control
ChR2_rVal,ChR2_p_raw,ChR2_p_adj, ChR2_alpha = corrMatrix(ChR2_raw_data)
Control_rVal,Control_p_raw,Control_p_adj, Control_alpha = corrMatrix(Control_raw_data)

#After getting the rVal and adjust p-values, then run the function to check for significance and generate the corr matrices with all non-zero values
ChR2_threshold_matrix,ChR2_pandas_matrix = significanceCheck(ChR2_p_adj, ChR2_rVal, alpha=ChR2_alpha, threshold = 0.001,
                                                             names=ChR2_nodes, plot=True, include_Negs=True, Anatomy=ROIs)
Control_threshold_matrix,Control_pandas_matrix = significanceCheck(Control_p_adj,Control_rVal,alpha = Control_alpha, threshold = 0.001,
                                                                   names=Control_nodes,plot = True,include_Negs = True, Anatomy=ROIs)

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph, ChR2_pos = networx(ChR2_threshold_matrix,ChR2_nodes)
Control_graph, Control_pos = networx(Control_threshold_matrix,Control_nodes)

#Run some hierarchical clustering
ChR2_hc_clusters,ChR2_hc_df= hierarch_clust(ChR2_graph,ChR2_nodes,ROIs.values(),plot = False)
Control_hc_clusters,Control_hc_df= hierarch_clust(Control_graph,Control_nodes,ROIs.values(),plot = False)

#Generate a quick plot for both groups showing cuts along the HC dendrogram
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(ChR2_hc_clusters['Distance'],ChR2_hc_clusters['Number of Clusters'],c='b',label='ChR2')
ax1.scatter(Control_hc_clusters['Distance'],Control_hc_clusters['Number of Clusters'],c='g',label='Control')
plt.xlabel("Cut Distance Along Dendrogram")
plt.ylabel("Number of Clusters")
plt.legend(loc='upper right')

#Run some markov clustering
ChR2_markov_df,ChR2_markov_clusters = markov(ChR2_graph,plot = False)
Control_markov_df,Control_markov_clusters = markov(Control_graph,plot = False)

