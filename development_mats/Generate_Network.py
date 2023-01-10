#Unpickle the ROI dictionary
import pickle as pkl
from random import shuffle
with open('/Users/kedorst/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle','rb') as f:
    ROIs = pkl.load(f)
Allen_Groups = list(ROIs.values())

#Then get that data
ChR2_raw_data, ChR2_nodes = loadData('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/ChR2_Small_Box.csv')
Control_raw_data, Control_nodes = loadData('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/Control_Small_Box.csv')

#Function to compare densities of the raw data only within areas across two conditions using a ttest
df_stats = comp_conds(ChR2_nodes, ChR2_raw_data, Control_raw_data)

#Get the correlation and adjusted p values for data_ChR2 and data_Control
ChR2_rVal, ChR2_p_Val = corrMatrix(ChR2_raw_data, corr_type="Spearman")
Control_rVal, Control_p_Val = corrMatrix(Control_raw_data, corr_type="Spearman")

#Can use the significanceCheck function to still get the raw correlation matrices
ChR2_threshold_matrix, ChR2_pandas_matrix = significanceCheck(ChR2_p_Val, ChR2_rVal, 0.05, names=ChR2_nodes, plot=True, Anatomy=ROIs)
Control_threshold_matrix, Control_pandas_matrix = significanceCheck(Control_p_Val, Control_rVal, 0.05, names=Control_nodes, plot=True, Anatomy=ROIs)

#A threshold array to look at the top 20% magnitude edges
ChR2_percent_array = percentile(ChR2_threshold_matrix, .25)
Control_percent_array = percentile(Control_threshold_matrix, .25)

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph = networx(ChR2_percent_array, ChR2_nodes)
Control_graph = networx(Control_percent_array, Control_nodes)

#Get the node attributes of the generated graphs
ChR2_node_attrs = grab_node_attributes(ChR2_graph, use_distance=False, compress_to_df=True)
Control_node_attrs = grab_node_attributes(Control_graph, use_distance=False, compress_to_df=True)

#Run some markov clustering
ChR2_markov_df, ChR2_markov_clusters, ChR2_markov_vector = markov(ChR2_graph, ChR2_nodes)
Control_markov_df, Control_markov_clusters, Control_markov_vector = markov(Control_graph, Control_nodes)

#Run some Louvain clustering
ChR2_lou_comm, ChR2_lou_max_mod, ChR2_lou_mod_mean, ChR2_lou_vector = louvain(ChR2_graph, ChR2_nodes, 1000)
shuffle(ChR2_lou_comm)
Control_lou_comm, Control_lou_max_mod, Control_lou_mod_mean, Control_lou_vector = louvain(Control_graph, Control_nodes, 1000)
shuffle(Control_lou_comm)

#Make visuals of the graphs using the set of plot_utils functions
#Set some parameters for spacing the nodes and getting colors according to Allen ROIs
allen_color_list = get_allen_colors('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/ROIs.csv')
sunflower = sunflower_theta(147)
sunflower_radius = sunflower_r(147)
point_cloud = get_point_cloud(k=0)


#the positions for plotting in a graph depend on the clustering algorithm performed
ChR2_pos_dict = get_position_data(ChR2_lou_comm, ChR2_nodes, shape=False)
Control_pos_dict = get_position_data(Control_lou_comm, Control_nodes, shape=False)


#Graph the things
graph_network(ChR2_graph, allen_color_list, pos_dict=ChR2_pos_dict)
graph_network(Control_graph, allen_color_list, pos_dict=Control_pos_dict)

#Based on the clustering algorithm, find the community that houses the DG
ChR2_DG_subgraph = DG_subgraph(ChR2_lou_comm, ChR2_nodes, ChR2_graph, ChR2_pos_dict, allen_color_list)
Control_DG_subgraph = DG_subgraph(Control_lou_comm, Control_nodes, Control_graph, Control_pos_dict, allen_color_list)

#Can also get the subgraph attributes
ChR2_DG_sub_attrs = grab_node_attributes(ChR2_DG_subgraph, use_distance=False, compress_to_df=True)
Control_DG_sub_attrs = grab_node_attributes(Control_DG_subgraph, use_distance=False, compress_to_df=True)

#For clustering, you can calculate the PC and WMDz for each node
ChR2_lou_WMDz_PC_df = cluster_attributes(ChR2_graph, ChR2_nodes, ChR2_lou_vector, make_df=True)
Control_lou_WMDz_PC_df = cluster_attributes(Control_graph, Control_nodes, Control_lou_vector, make_df=True)


#Find some hub ROIs based on metrics from Van den Huevel(2010)
ChR2_whole_graph_Results, ChR2_whole_graph_Hubs = findMyHubs(ChR2_node_attrs)
Control_whole_graph_Results, Control_whole_graph_Hubs = findMyHubs(Control_node_attrs)

ChR2_DGsub_Results, ChR2_DGsub_Hubs = findMyHubs(ChR2_DG_sub_attrs)
Control_DGsub_Results, Control_DGsub_Hubs = findMyHubs(Control_DG_sub_attrs)

#Optional, combine all of the attributes into one final dataframe
ChR2_final_df = combine_node_attrs(ChR2_whole_graph_Results, ChR2_lou_WMDz_PC_df, Allen_Groups)
Control_final_df = combine_node_attrs(Control_whole_graph_Results, Control_lou_WMDz_PC_df, Allen_Groups)

#If you wish to export all of your data to .csv files, run the node_attrs_to_csv function
node_attrs_to_csv(ChR2_final_df, '/Users/kedorst/Desktop/ChR2_Large_Box', 'ChR2_nodes_Large_Box')
node_attrs_to_csv(Control_final_df, '/Users/kedorst/Desktop/Control_Large_Box', 'Control_nodes_Large_Box')
