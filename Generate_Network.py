#Start by loading the functions from the referenced NetworkFunctions.py file
run NetworkFunctions.py

#Start by getting the data
ChR2_data = '/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/ChR2/ChR2_Large_Network.csv'
raw_data_ChR2,nodes_ChR2 = loadData(ChR2_data)

Control_data = '/Users/kaitlyndorst/Desktop/Data_Analyses/Networks/Network Wrangling/Control/Control_Large_Network.csv'
raw_data_Control,nodes_Control = loadData(Control_data)

#Get the correlation and adjusted p values for data_ChR2 and data_Control
rVal_ChR2,p_ChR2,p_adj_ChR2, ChR2_alpha = corrMatrix(raw_data_ChR2)
rVal_Control,p_Control,p_adj_Control, Control_alpha = corrMatrix(raw_data_Control)

#After getting the rVal and adjust p-values, then run the function to check for significance and generate the corr matrices
ChR2_threshold_matrix, ChR2_pandas_matrix = significanceCheck(p_adj_ChR2, rVal_ChR2, alpha = ChR2_alpha, threshold = 0.001, names=nodes_ChR2, plot = True, include_Negs=True)
plt.savefig("ChR2_corr_matrix_thresh0_ward.jpeg",optimize = True)

Control_threshold_matrix, Control_pandas_matrix = significanceCheck(p_adj_Control,rVal_Control,alpha = Control_alpha, threshold = 0.001, names=nodes_Control,plot = True,include_Negs = True)
plt.savefig("Control_corr_matrix_thresh0_ward.jpeg",optimize = True)

#Run some hierarchical clustering

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph, ChR2_pos = networx(ChR2_threshold_matrix,nodes_ChR2)
Control_graph, Control_pos = networx(Control_threshold_matrix,nodes_Control)

#Run some markov clustering
