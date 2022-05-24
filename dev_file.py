import network_analysis
import network_analysis as net_a
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx



if __name__ == "__main__":
    # THIS IS FOR TEST PURPOSES ONLY
    file = '/home/ryansenne/PycharmProjects/Networks2/csv_files/ChR2_Large_Network.csv'
    file2 = '/home/ryansenne/PycharmProjects/Networks2/csv_files/Control_Large_Network.csv'
    # allen_groups = pd.read_csv('/Users/ryansenne/PycharmProjects/networkx/csv_files/ROIs.csv')
    test_data, test_nodes = net_a.loadData(file)
    rvals, p= net_a.corrMatrix(test_data, 'Spearman')
    threshold_matrix = net_a.significanceCheck(p, rvals, 1, names=test_nodes, include_Negs=False)
    per = net_a.percentile(threshold_matrix, 0.2)
    G, pos = net_a.networx(per, test_nodes)

    test_data1, test_nodes1 = net_a.loadData(file2)
    rvals1, p1= net_a.corrMatrix(test_data1, corr_type='Spearman')
    threshold_matrix1 = net_a.significanceCheck(p1, rvals1, 1, names=test_nodes1, include_Negs=False)
    per1 = net_a.percentile(threshold_matrix1, 0.2)
    G1, pos1 = net_a.networx(per1, test_nodes1)

    df, mark_clust, _ = net_a.markov(G, test_nodes)
    df1, mark_clust1, _1 = net_a.markov(G1, test_nodes1)

    pos_dict = net_a.get_position_data(mark_clust, test_nodes)
    my_allen_colors = network_analysis.get_allen_colors('ROIs.csv')
    net_a.graph_network(G, my_allen_colors, pos_dict)

    pos_dict1 = net_a.get_position_data(mark_clust, test_nodes)
    my_allen_colors1 = network_analysis.get_allen_colors('ROIs.csv')
    net_a.graph_network(G1, my_allen_colors1, pos_dict1)
    # # color_list = net_a.grab_color_attributes(mark_clust, test_nodes)
    # # my_del = net_a.in_silico_deletion(G, plot=True)
    # my_list = net_a.get_ordered_degree_list(G)
    # pe, m = net_a.threshold_simulation(threshold_matrix, 0.05, 0.7, 10)
    # pe1, m1 = net_a.threshold_simulation(threshold_matrix1, 0.05, 0.7, 10)
    # clust, result = net_a.hierarch_clust(threshold_matrix, test_nodes, allen_groups['Allen Group Name'])
    # plt.plot(pe, m)
    # plt.plot(pe1, m1)
