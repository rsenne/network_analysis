#Unpickle the ROI dictionary
import pickle as pkl
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

#Get the node attributes of the generated graphs
ChR2_node_degree,ChR2_node_eig,ChR2_node_between,ChR2_node_close = grab_node_attributes(ChR2_graph)
Control_node_degree,Control_node_eig,Control_node_between,Control_node_close = grab_node_attributes(Control_graph)

#If you wish to export all of your data to .csv files, run the node_attrs_to_csv function
ChR2_str_nodes = nameof(ChR2_nodes)
ChR2_nodes_df= node_attrs_to_csv(ChR2_nodes,ChR2_node_degree,ChR2_node_between,ChR2_node_eig,ChR2_node_close,'/Users/kaitlyndorst/Desktop/Network_csv',ChR2_str_nodes)

Control_str_nodes = nameof(Control_nodes)
Control_nodes_df = node_attrs_to_csv(Control_nodes,Control_node_degree,Control_node_between,Control_node_eig,Control_node_close,'/Users/kaitlyndorst/Desktop/Network_csv',Control_str_nodes)

#Run some hierarchical clustering
ChR2_hc_cuts_df,ChR2_hc_assigns,ChR2_hc_clusters= hierarch_clust(ChR2_graph,ChR2_nodes,ROIs.values(),plot = False)
Control_hc_cuts_df,Control_hc_assigns,Control_hc_clusters= hierarch_clust(Control_graph,Control_nodes,ROIs.values(),plot = False)

#Generate a quick plot for both groups showing cuts along the HC dendrogram
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(ChR2_hc_cuts_df['Distance Cut'],ChR2_hc_cuts_df['Number of Clusters'],c='b',label='ChR2')
ax1.scatter(Control_hc_cuts_df['Distance Cut'],Control_hc_cuts_df['Number of Clusters'],c='g',label='Control')
plt.xlabel("Cut Distance Along Dendrogram")
plt.ylabel("Number of Clusters")
plt.legend(loc='upper right')

#Run some markov clustering
ChR2_markov_df,ChR2_markov_clusters = markov(ChR2_graph,plot = False)
Control_markov_df,Control_markov_clusters = markov(Control_graph,plot = False)

#Run louvain clustering
ChR2_lou_clust = louvain(ChR2_graph,ChR2_nodes)
Control_lou_clust = louvain(Control_graph,Control_nodes)
