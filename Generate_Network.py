#Start by loading the functions from the referenced NetworkFunctions.py file
run NetworkFunctions.py

#Need some ROI info again
ROIs = pd.read_csv("/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/ROIs.csv")
ROIs = ROIs.loc[:,["Abbreviation","Allen Group Name"]].sort_values("Allen Group Name").reset_index(drop = True)

#Then get that data
ChR2_data = '/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/ChR2/ChR2_Large_Network.csv'
ChR2_raw_data,ChR2_nodes = loadData(ChR2_data)

Control_data = '/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/Control/Control_Large_Network.csv'
Control_raw_data,Control_nodes = loadData(Control_data)

#Get the correlation and adjusted p values for data_ChR2 and data_Control
ChR2_rVal,ChR2_p_raw,ChR2_p_adj, ChR2_alpha = corrMatrix(ChR2_raw_data)
Control_rVal,Control_p_raw,Control_p_adj, Control_alpha = corrMatrix(Control_raw_data)

#After getting the rVal and adjust p-values, then run the function to check for significance and generate the corr matrices with all non-zero values
ChR2_threshold_matrix = significanceCheck(ChR2_p_adj, ChR2_rVal, alpha = ChR2_alpha, threshold = 0.001, names=ChR2_nodes, plot = False, include_Negs=True)
Control_threshold_matrix = significanceCheck(Control_p_adj,Control_rVal,alpha = Control_alpha, threshold = 0.001, names=Control_nodes,plot = False,include_Negs = True)

#Run some hierarchical clustering
ChR2_hc_clusters,ChR2_hc_df,ChR2_components = hierarch_clust(ChR2_threshold_matrix,ChR2_nodes,ROIs['Allen Group Name'],plot = False)
Control_hc_clusters,Control_hc_df,Control_components = hierarch_clust(Control_threshold_matrix,Control_nodes,ROIs['Allen Group Name'],plot = False)

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph, ChR2_pos = networx(ChR2_threshold_matrix,ChR2_nodes)
Control_graph, Control_pos = networx(Control_threshold_matrix,Control_nodes)

#Run some markov clustering
ChR2_markov_df,ChR2_markov_results,ChR2_markov_clusters = markov(ChR2_graph,plot = False)
plt.savefig("ChR2_markov_graph.jpeg",optimize = True)

Control_markov_df,Control_markov_results,Control_markov_clusters = markov(Control_graph,plot = False)
plt.savefig("Control_markov_graph.jpeg",optimize = True)

