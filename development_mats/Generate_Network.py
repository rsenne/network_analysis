#Unpickle the ROI dictionary
import pickle as pkl
from scipy import stats
with open('/Users/kedorst/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle','rb') as f:
    ROIs = pkl.load(f)
Allen_Groups = list(ROIs.values())

#Then get that data
ChR2_raw_data,ChR2_nodes = loadData('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/ChR2_Small_Box.csv')
Control_raw_data,Control_nodes = loadData('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/Control_Small_Box.csv')

#Function to compare densities of the raw data using a Kruskal-Wallis H test
df_stats = comp_conds(ChR2_nodes,ChR2_raw_data,Control_raw_data)

#Get the correlation and adjusted p values for data_ChR2 and data_Control
ChR2_rVal,ChR2_p_Val = corrMatrix(ChR2_raw_data, corr_type="Spearman")
Control_rVal,Control_p_Val = corrMatrix(Control_raw_data, corr_type="Spearman")

#A threshold array to look at the top 20% magnitude edges
ChR2_percent_array = percentile(ChR2_rVal,.20)
Control_percent_array = percentile(Control_rVal,.20)

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph = networx(ChR2_percent_array,ChR2_nodes)
Control_graph = networx(Control_percent_array,Control_nodes)

#Get the node attributes of the generated graphs
ChR2_node_attrs = grab_node_attributes(ChR2_graph,use_distance=False,compress_to_df=True)
Control_node_attrs = grab_node_attributes(Control_graph,use_distance=False,compress_to_df=True)


#Run some markov clustering
ChR2_markov_df,ChR2_markov_clusters,ChR2_markov_vector = markov(ChR2_graph,ChR2_nodes)
Control_markov_df,Control_markov_clusters,Control_markov_vector = markov(Control_graph,Control_nodes)


#Make visuals of the graphs using the set of plot_utils functions
#Set some parameters for spacing the nodes and getting colors according to Allen ROIs
allen_color_list = get_allen_colors('/Users/kedorst/Documents/GitHub/network_analysis/csv_files/ROIs.csv')
sunflower = sunflower_theta(147)
sunflower_radius = sunflower_r(147)
point_cloud = get_point_cloud(k=0)


#the positions for plotting in a graph depend on the clustering algorithm performed
ChR2_pos_dict = get_position_data(ChR2_markov_clusters,ChR2_nodes)
Control_pos_dict = get_position_data(Control_markov_clusters,Control_nodes)


#Graph the things
graph_network(ChR2_graph,allen_color_list,pos_dict=ChR2_pos_dict)
graph_network(Control_graph,allen_color_list,pos_dict=Control_pos_dict)


#For clustering, you can calculate the PC and WMDz for each node
ChR2_mc_WMDz_PC_df = cluster_attributes(ChR2_graph,ChR2_nodes,ChR2_markov_vector,make_df=True)
Control_mc_WMDz_PC_df = cluster_attributes(Control_graph,Control_nodes,Control_markov_vector,make_df=True)


#Find some hub ROIs based on metrics from Van den Huevel(2010)
ChR2_Results,ChR2_Hubs = findMyHubs(ChR2_node_attrs)
Control_Results,Control_Hubs = findMyHubs(Control_node_attrs)


#Optional, combine all of the attributes into one final dataframe
ChR2_final_df = combine_node_attrs(ChR2_Results,ChR2_mc_WMDz_PC_df,Allen_Groups,ChR2_delta_glob_eff)
Control_final_df = combine_node_attrs(Control_Results,Control_mc_WMDz_PC_df,Allen_Groups,Control_delta_glob_eff)

#If you wish to export all of your data to .csv files, run the node_attrs_to_csv function
node_attrs_to_csv(ChR2_final_df,'/Users/kaitlyndorst/Desktop/ChR2_Small_Box','ChR2_nodes_Small_Box')
node_attrs_to_csv(Control_final_df,'/Users/kaitlyndorst/Desktop/Control_Small_Box','Control_nodes_Small_Box')
