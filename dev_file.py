import network_analysis
import network_analysis as net_a
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import multiprocessing as mp


if __name__ == "__main__":
    # THIS IS FOR TEST PURPOSES ONLY
    control_of = r'C:\Users\Ryan Senne\PycharmProjects\network_analysis\csv_files\ChR2_Small_Box.csv'
    ChR2_of = r'C:\Users\Ryan Senne\PycharmProjects\network_analysis\csv_files\Control_Small_Box.csv'
    # allen_groups = pd.read_csv('/Users/ryansenne/PycharmProjects/networkx/csv_files/ROIs.csv')
    test_data, test_nodes = net_a.loadData(control_of)
    rvals, p = net_a.corrMatrix(test_data, corr_type='Spearman', z_trans=False)
    threshold_matrix = net_a.significanceCheck(p, rvals, 1, names=test_nodes, include_Negs=True)
    per = net_a.percentile(threshold_matrix, 0.2)
    G = net_a.networx(per, test_nodes)

    test_data1, test_nodes1 = net_a.loadData(ChR2_of)
    rvals1, p1 = net_a.corrMatrix(test_data1, corr_type='Spearman', z_trans=False)
    threshold_matrix1 = net_a.significanceCheck(p1, rvals1, 1, names=test_nodes1, include_Negs=True)
    per1 = net_a.percentile(threshold_matrix1, 0.2)
    G1 = net_a.networx(per1, test_nodes1)

    df, mark_clust, p8, clust_dit = net_a.markov(G, test_nodes)
    df1, mark_clust1, _, clust_dict_1 = net_a.markov(G1, test_nodes1)
    #
    pos_dict = net_a.get_position_data(mark_clust, test_nodes, shape=False)
    my_allen_colors = network_analysis.get_allen_colors(r'C:\Users\Ryan Senne\PycharmProjects\network_analysis\csv_files\ROIs.csv')
    # r1 = net_a.graph_network(G, my_allen_colors, pos_dict)

    pos_dict1 = net_a.get_position_data(mark_clust1, test_nodes1, shape=False)
    my_allen_colors = network_analysis.get_allen_colors(
        r'C:\Users\Ryan Senne\PycharmProjects\network_analysis\csv_files\ROIs.csv')
    # r = net_a.graph_network(G1, my_allen_colors, pos_dict1)

    chr2_hub, top_chr2 = net_a.findMyHubs(net_a.grab_node_attributes(G1, compress_to_df=True))
    control_hub, top_control = net_a.findMyHubs(net_a.grab_node_attributes(G, compress_to_df=True))

    control_of_GE = nx.global_efficiency(G)
    ChR2_of_GE = nx.global_efficiency(G1)

    control_of_LE = nx.local_efficiency(G)
    ChR2_of_LE = nx.local_efficiency(G1)



    # control_of_Sigma = nx.sigma(G)
    # ChR2_of_Sigma = nx.sigma(G1)
