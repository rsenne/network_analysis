#Unpickle the ROI dictionary
import pickle as pkl
from scipy import stats
with open('/Users/kaitlyndorst/Documents/GitHub/networkx/Allen_Areas_dict.pickle','rb') as f:
    ROIs = pkl.load(f)
Allen_Groups = list(ROIs.values())

#Then get that data
ChR2_raw_data,ChR2_nodes = loadData('/Users/kaitlyndorst/Desktop/ChR2_Small_Box/ChR2_Small_Box.csv')
Control_raw_data,Control_nodes = loadData('/Users/kaitlyndorst/Desktop/Control_Small_Box/Control_Small_Box.csv')

#Function to compare densities of the raw data using a Kruskal-Wallis H test
df_stats = comp_conds(ChR2_nodes,ChR2_raw_data,Control_raw_data)

#Get the correlation and adjusted p values for data_ChR2 and data_Control
ChR2_rVal,ChR2_p_raw,ChR2_p_adj,ChR2_alpha = corrMatrix(ChR2_raw_data)
Control_rVal,Control_p_raw,Control_p_adj,Control_alpha = corrMatrix(Control_raw_data)

#After getting the rVal and adjust p-values, then run the function to check for significance and generate the corr matrices with all non-zero values
ChR2_threshold_matrix,ChR2_pandas = significanceCheck(ChR2_p_adj, ChR2_rVal, alpha=ChR2_alpha, threshold=0.0,
                                                             names=ChR2_nodes, plot=True, include_Negs=True, Anatomy=ROIs)
Control_threshold_matrix,Control_pandas = significanceCheck(Control_p_adj,Control_rVal,alpha=Control_alpha, threshold=0.0,
                                                                   names=Control_nodes, plot=True, include_Negs=True, Anatomy=ROIs)

#A threshold array to look at the top 20% magnitude edges
ChR2_percent_array = percentile(ChR2_threshold_matrix,.10)
Control_percent_array = percentile(Control_threshold_matrix,.10)

#Using the generated correlation matrices, build a graph using networkx
ChR2_graph,ChR2_pos = networx(ChR2_percent_array,ChR2_nodes)
Control_graph,Control_pos = networx(Control_percent_array,Control_nodes)

#Run some hierarchical clustering
ChR2_hc_cuts_df,ChR2_hc_mods,ChR2_hc_assigns,ChR2_hc_clusters= hierarch_clust(ChR2_graph,ChR2_nodes,ROIs.values(),plot = False)
Control_hc_cuts_df,Control_hc_mods,Control_hc_assigns,Control_hc_clusters= hierarch_clust(Control_graph,Control_nodes,ROIs.values(),plot = False)

#Generate a quick plot for both groups showing cuts along the HC dendrogram
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(ChR2_hc_cuts_df['Distance Cut'],ChR2_hc_cuts_df['Number of Clusters'],c='b',label='ChR2')
ax1.scatter(Control_hc_cuts_df['Distance Cut'],Control_hc_cuts_df['Number of Clusters'],c='g',label='Control')
plt.xlabel("Cut Distance Along Dendrogram")
plt.ylabel("Number of Clusters")
plt.legend(loc='upper right')

#Run some markov clustering
ChR2_markov_df,ChR2_markov_clusters,ChR2_markov_vector = markov(ChR2_graph,ChR2_nodes)
Control_markov_df,Control_markov_clusters,Control_markov_vector = markov(Control_graph,Control_nodes)

#With markov clustering, you can pass the clustering algorithm into a threshold function to examine how modularity is affected by edges
ChR2_mc_percentiles,ChR2_mc_modularity = threshold_simulation(ChR2_threshold_matrix,.05,.20,5)
Control_mc_percentiles,Control_mc_modularity = threshold_simulation(Control_threshold_matrix,.05,.20,5)

#Run louvain clustering
ChR2_lou_clust = louvain(ChR2_graph,ChR2_nodes)
Control_lou_clust = louvain(Control_graph,Control_nodes)

#Make visuals of the graphs using the set of plot_utils functions
#Set some parameters for spacing the nodes and getting colors according to Allen ROIs
allen_color_list = get_allen_colors('/Users/kaitlyndorst/Documents/GitHub/networkx/csv_files/ROIs.csv')
sunflower = sunflower_theta(155)
sunflower_radius = sunflower_r(155,c=1.5)
point_cloud = get_point_cloud(k=0)

#the positions for plotting in a graph depend on the clustering algorithm performed
ChR2_pos_dict = get_position_data(ChR2_markov_clusters,ChR2_nodes)
Control_pos_dict = get_position_data(Control_markov_clusters,Control_nodes)

#Graph the things
graph_network(ChR2_graph,allen_color_list,pos_dict=ChR2_pos_dict)
graph_network(Control_graph,allen_color_list,pos_dict=Control_pos_dict)

#Get the node attributes of the generated graphs
ChR2_node_attrs = grab_node_attributes(ChR2_graph,use_distance=False,compress_to_df=True)
Control_node_attrs = grab_node_attributes(Control_graph,use_distance=False,compress_to_df=True)

#For clustering, you can calculate the PC and WMDz for each node
ChR2_mc_WMDz_PC_df = cluster_attributes(ChR2_graph,ChR2_nodes,ChR2_markov_vector,make_df=True)
Control_mc_WMDz_PC_df = cluster_attributes(Control_graph,Control_nodes,Control_markov_vector,make_df=True)

#Find some hub ROIs based on metrics from Van den Huevel(2010)
ChR2_Results,ChR2_Hubs = findMyHubs(ChR2_node_attrs)
Control_Results,Control_Hubs = findMyHubs(Control_node_attrs)

#Do some quick in silico deletions and measure global efficiency
ChR2_delta_glob_eff = in_silico_deletion(ChR2_graph,plot=True)
Control_delta_glob_eff = in_silico_deletion(Control_graph,plot=True)

#Optional, combine all of the attributes into one final dataframe
ChR2_final_df = combine_node_attrs(ChR2_Results,ChR2_mc_WMDz_PC_df,Allen_Groups,ChR2_delta_glob_eff)
Control_final_df = combine_node_attrs(Control_Results,Control_mc_WMDz_PC_df,Allen_Groups,Control_delta_glob_eff)

#If you wish to export all of your data to .csv files, run the node_attrs_to_csv function
node_attrs_to_csv(ChR2_final_df,'/Users/kaitlyndorst/Desktop/ChR2_Small_Box','ChR2_nodes_Small_Box')
node_attrs_to_csv(Control_final_df,'/Users/kaitlyndorst/Desktop/Control_Small_Box','Control_nodes_Small_Box')

#Take a holistic approach and look at all of the edges in the network instead of the strongest
whole_ChR2_graph,whole_ChR2_pos = networx(ChR2_threshold_matrix,ChR2_nodes)
graph_network(whole_ChR2_graph,allen_color_list,whole_ChR2_pos)

whole_ChR2_node_attrs = grab_node_attributes(whole_ChR2_graph,use_distance=False,compress_to_df=True)
whole_ChR2_degree_list = get_ordered_degree_list(whole_ChR2_graph)
whole_ChR2_Results, whole_ChR2_Hubs = findMyHubs(whole_ChR2_node_attrs)

#Disruption Propagation Example
ChR2_VTA_FinalMat = disruptPropagate(whole_ChR2_graph,'VTA')
ChR2_VTA_FinalMat = ChR2_VTA_FinalMat.to_numpy()
VTA_del_graph,VTA_pos_dict = networx(ChR2_VTA_FinalMat,ChR2_nodes)
graph_network(VTA_del_graph,allen_color_list,pos_dict=whole_ChR2_pos)