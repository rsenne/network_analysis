import network_analysis as net_a
import pandas as pd
import networkx as nx

# THIS IS FOR TEST PURPOSES ONLY
file = '/home/ryansenne/PycharmProjects/Networks/csv_files/ChR2_Large_Network.csv'
file2 = '/home/ryansenne/PycharmProjects/Networks/csv_files/Control_Large_Network.csv'
allen_groups = pd.read_csv('/home/ryansenne/PycharmProjects/Networks/csv_files/ROIs.csv')

test_data, test_nodes = net_a.loadData(file)
rvals, p, p_adj, rej = net_a.corrMatrix(test_data)
threshold_matrix = net_a.significanceCheck(p_adj, rvals, 0.001, names=test_nodes, include_Negs=False)
G, pos = net_a.networx(threshold_matrix, test_nodes)
df, mark_clust = net_a.markov(G)
color_list = net_a.grab_color_attributes(mark_clust, test_nodes)
pos_dict = net_a.get_position_data(mark_clust, test_nodes)
net_a.graph_network(G, list(color_list.values()), pos_dict)
# my_del = net_a.in_silico_deletion(G, plot=True)
my_list = net_a.get_ordered_degree_list(G)
# clust, result = net_a.hierarch_clust(threshold_matrix, test_nodes, allen_groups['Allen Group Name'])