"""
this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by Drs. Cesar Coelho, Anne Wheeler, and Gisella Vetere.
we thank them for their kind support throughout this process
"""
# author:ryan senne/ramirez group

# import necessary libraries

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.stats
import seaborn as sns
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.patches as mpatches
# from bct.algorithms import centrality
from scipy.spatial.distance import cdist
# import network_analysis.algorithms as algorithms


# simple function for loading our csv file
def loadData(data):
    """
    :param data:
    :return:
    """
    data = pd.read_csv(data)
    data = data.apply(lambda x: x.fillna(x.mean()), axis=0)
    node_names = data.columns.to_list()
    node_number = list(item for item in range(0, len(node_names)))
    nodes = {node_number[i]: node_names[i] for i in range(len(node_number))}
    return data, nodes


def comp_conds(nodes, data1, data2):
    nodes = list(nodes.values())
    H_stat = []
    p_val = []
    for node in nodes:
        array1 = np.array(data1[node])
        array2 = np.array(data2[node])


        H, p = scipy.stats.kruskal(array1, array2)

        H_stat.append(H)
        p_val.append(p)

    d = {"H statistic": H_stat, "p value": p_val}
    df_stats = pd.DataFrame(d, columns=["H statistic", "p value"], index=nodes)
    df_stats["significant?"] = np.where(df_stats["p value"] < 0.05, True, False)

    return df_stats


# correlate our c-Fos counts between brain regions, df for data
# type for correlation coefficient i.e. "pearson"
def corrMatrix(data, corr_type='Pearson', z_trans=True):
    if corr_type == 'Pearson':
        rVal = np.corrcoef(data, rowvar=False)  # calculate pearson coefficients
        rVal[np.isnan(rVal)] = 0  # Will make all NaN values into zero
        rf = rVal[np.triu_indices(rVal.shape[0], 1)]  # upper triangular matrix of data to shuffle for p-value calc
        df = data.shape[1] - 2  # calculate degrees of freedom
        ts = rf * rf * (df / (1 - rf * rf))  # calculate t's
        pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))  # calculate p's from beta incomplete function
        # generate p-value matrix
        p = np.zeros(shape=rVal.shape)
        p[np.triu_indices(p.shape[0], 1)] = pf
        p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]
        p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
        # Multiple comparison of p values using Bonferroni correction
        rejected, p_adjusted, _, alpha_corrected = multipletests(p, alpha=0.05, method='bonferroni', is_sorted=True)
        np.fill_diagonal(rVal, 0)  # set main diagonal zero to avoid errors
        if z_trans:
            return np.arctanh(rVal), p, p_adjusted, alpha_corrected
        else:
            return rVal, p, p_adjusted, alpha_corrected
    elif corr_type == 'Spearman':
        rVal, pVal = scipy.stats.spearmanr(data, axis=0)
        return rVal, pVal


# Function that will threshold R-vals on the top percentile
def percentile(array, p):
    num_obs = int(np.size(array, 0) ** 2 * p)
    crit_value = -np.sort(-array.flatten())[num_obs - 1]
    percent_arr = np.where(array < crit_value, 0, array)
    return percent_arr


# Will generate a euclidean distance matrix from the raw data
def euclMatrix(data):
    data = data.T
    eucl_matrix = cdist(data, data, metric='euclidean')
    return eucl_matrix


# using this function we will threshold based off of p-values previously calculated
def significanceCheck(p_adjusted, corr, alpha, threshold=0.0, names=None, plot=False, include_Negs=True, Anatomy=None):
    p_adjusted = np.where((p_adjusted >= alpha), 0, p_adjusted)  # if not significant --> zero
    np.fill_diagonal(p_adjusted, 0)
    p_adjusted = np.where((p_adjusted != 0), 1, p_adjusted)  # if significant --> one
    if not include_Negs:
        p_adjusted = np.where((corr < 0), 0, p_adjusted)
    threshold_matrix = np.multiply(corr, p_adjusted)  # remove any insignificant correlations
    # remove correlations below threshold
    threshold_matrix = np.where((abs(threshold_matrix) < threshold), 0, threshold_matrix)
    # create a heatmap of correlations if wanted
    if plot:
        if names:
            pandas_matrix = pd.DataFrame(threshold_matrix, index=list(names.values()), columns=list(names.values()))
            if Anatomy:
                # Create a sorted dictionary from the unpickled ROIs dictionary
                sorted_dict = dict(sorted(Anatomy.items(), key=lambda item: item[1]))
                list_for_sort = list(sorted_dict.keys())

                # Reorganize the pandas_matrix to reflect the order of the Allen ROIs
                allen_pandas = pandas_matrix[list_for_sort].loc[list_for_sort]

                # Create color assignments for Allen ROIs
                num_allens = list(sorted_dict.values())
                allens_unique = np.unique(num_allens)
                color_list = [color for color in sns.color_palette('Set3', n_colors=len(allens_unique))]
                color_ref = dict(zip(map(str, allens_unique), color_list))
                allen_colors = pd.Series(num_allens, index=allen_pandas.columns).map(color_ref)

                # Create a legend for the Allen ROIs
                cerebellum = mpatches.Patch(color=color_list[0], label='Cerebellum')
                cort_plate = mpatches.Patch(color=color_list[1], label='Cortical Plate')
                cort_subplate = mpatches.Patch(color=color_list[2], label='Cortical Subplate')
                hypothalamus = mpatches.Patch(color=color_list[3], label='Hypothalamus')
                medulla = mpatches.Patch(color=color_list[4], label='Medulla')
                midbrain = mpatches.Patch(color=color_list[5], label='Midbrain')
                pallidum = mpatches.Patch(color=color_list[6], label='Pallidum')
                pons = mpatches.Patch(color=color_list[7], label='Pons')
                striatum = mpatches.Patch(color=color_list[8], label='Striatum')
                thalamus = mpatches.Patch(color=color_list[9], label='Thalamus')

                # Plot the newly generated Allen ROI correlation maitrx
                plt.figure()
                sns.clustermap(allen_pandas, cmap='vlag', row_colors=allen_colors, col_colors=allen_colors,
                               row_cluster=False, col_cluster=False, xticklabels=False, yticklabels=False,
                               figsize=(10, 10), cbar_pos=(0.1, 0.15, .02, .4), cbar_kws={'label': 'arctan(R)'})
                plt.legend(
                    handles=[cerebellum, cort_plate, cort_subplate, hypothalamus, medulla, midbrain, pallidum, pons,
                             striatum, thalamus],
                    bbox_to_anchor=(5.0, 1.6))
        else:
            pandas_matrix = pd.DataFrame(threshold_matrix)
        sns.clustermap(pandas_matrix, cmap='vlag', method='ward', metric='euclidean', figsize=(10, 10),
                       cbar_pos=(.9, .9, .02, .10))
        return threshold_matrix, pandas_matrix
    else:
        return threshold_matrix


# we will create our undirected network graphs based on our matrices
def networx(corr_data, nodeLabel, drop_islands=False):
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)  # use the updated corr_data to make a graph
    if nodeLabel:
        graph = nx.relabel_nodes(graph, nodeLabel)
    graph = nx.relabel_nodes(graph, nodeLabel)
    if drop_islands:
        remove_node = [node for node, degree in graph.degree() if degree < 1]
        graph.remove_nodes_from(remove_node)
    return graph



def lazy_network_generator(data):
    df, nodes = loadData(data)
    rVal, p, p_adjusted, alpha_corrected = corrMatrix(df)
    threshold_matrix = significanceCheck(p_adjusted, rVal, alpha=0.001, names=nodes)
    G, pos = networx(threshold_matrix, nodes)
    return G


def shortest(G, threshold_matrix):
    # Function to calculate shortest path of each node
    short = nx.floyd_warshall_numpy(G, weight='weight')
    shortavg = np.mean(short, axis=0)
    keys = threshold_matrix.index.values.tolist()
    vals = shortavg.tolist()
    zip_iterator = zip(keys, vals)
    short_dictionary = dict(zip_iterator)

    return short_dictionary


def grab_node_attributes(graph, use_distance=False, compress_to_df=False):
    if use_distance:
        G_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in
                           graph.edges(data='weight')}  # creates a dict of calculted distance between all nodes
        nx.set_edge_attributes(graph, values=G_distance_dict, name='distance')
    deg = nx.degree_centrality(graph)
    between = nx.betweenness_centrality(graph)
    eig = nx.eigenvector_centrality(graph)
    close = nx.closeness_centrality(graph)
    clust = nx.clustering(graph)
    deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[0])}
    between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[0])}
    eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[0])}
    close_sort = {area: val for area, val in sorted(close.items(), key=lambda ele: ele[0])}
    clust_sort = {area: val for area, val in sorted(clust.items(), key=lambda ele: ele[0])}
    if compress_to_df:
        node_info = {
            'Degree': list(deg_sort.values()),
            'Betweenness': list(between_sort.values()),
            'Eigenvector_Centrality': list(eig_sort.values()),
            'Closeness': list(close_sort.values()),
            'Clustering_Coefficient': list(clust_sort.values()),
        }
        ROI_index = list(graph.nodes)
        node_attrs_df = pd.DataFrame(node_info, index=ROI_index, columns=node_info.keys())
        return node_attrs_df
    else:
        return deg_sort, between_sort, eig_sort, close_sort, clust_sort


def get_ordered_list(G, stat='Degree'):
    if stat == 'Degree':
        ordered = {k: v for k, v in sorted(dict(G.degree()).items(), key=lambda item: item[1])}
    elif stat == 'WeightedDegree':
        ordered = {k: v for k, v in sorted(dict(G.degree(weight='weight')).items(), key=lambda item: item[1])}
    elif stat == 'EigCentrality':
        ordered = {k: v for k, v in sorted(dict(nx.eigenvector_centrality(G)).items(), key=lambda item: item[1])}
    elif stat == 'Clustering':
        ordered = {k: v for k, v in sorted(dict(nx.clustering(G)).items(), key=lambda item: item[1])}
    elif stat == 'Betweenness':
        ordered = {k: v for k, v in sorted(dict(nx.betweenness_centrality(G)).items(), key=lambda item: item[1])}
    else:
        ordered = []
    return ordered


def cluster_attributes(graph, nodes, communities, make_df=False):
    adj_matrix = nx.to_numpy_array(graph)  # Will create an adjacency matrix from the graph as a np.ndarray
    node_ROIs = nodes.values()
    WMDz = centrality.module_degree_zscore(adj_matrix, communities, flag=0)  # calculate the WMDz
    PC = centrality.participation_coef(adj_matrix, communities, 'undirected')  # calculate the participation coefficient
    if make_df:
        d = {'WMDz': WMDz, "PC": PC}
        df = pd.DataFrame(d, columns=["WMDz", "PC"], index=node_ROIs)
        return df
    else:
        return WMDz, PC


def findMyHubs(node_attrs_df):
    Results = node_attrs_df
    # Van den Huevel(2010) - https://www.jneurosci.org/content/30/47/15915
    # used the top or bottom quartiles to determine the hubness of all nodes so here we calculate that.
    # For each significant measure an ROI has add one in the score column, a score >= 2 is considered a hub node.

    Results['Hub_Score'] = 0
    Results['Hub_Score'] = np.where((Results['Degree'] >= Results.Degree.quantile(0.80)), Results['Hub_Score'] + 1,
                                    Results['Hub_Score'])
    Results['Hub_Score'] = np.where((Results['Eigenvector_Centrality'] >= Results.Eigenvector_Centrality.quantile(.80)),
                                    Results['Hub_Score'] + 1,
                                    Results['Hub_Score'])
    Results['Hub_Score'] = np.where((Results['Betweenness'] >= Results.Betweenness.quantile(0.80)),
                                    Results['Hub_Score'] + 1,
                                    Results['Hub_Score'])
    Results['Hub_Score'] = np.where((Results['Clustering_Coefficient'] >= Results.Clustering_Coefficient.quantile(.80)),
                                    Results['Hub_Score'] + 1,
                                    Results['Hub_Score'])
    Results['Hub_Score'] = np.where((Results['Closeness'] >= Results.Closeness.quantile(.80)),       #here is the adding function for shortest path
                                    Results['Hub_Score'] + 1, Results['Hub_Score'])

    NonHubs = Results[
        (Results['Hub_Score'] < 2)].index  # create an index of rois with a score of less than 2 in hubness

    Hubs = Results.drop(NonHubs).sort_values('Hub_Score',
                                             ascending=False)  # create a new frame with only the important nodes/ take out rois in the prior index

    return Results, Hubs


def threshold_simulation(adj_mat, a, b, x, algo='markov'):
    percentiles = [i for i in np.linspace(a, b, x)]
    thresholded_arrays = [percentile(adj_mat, p) for p in percentiles]
    if algo == 'markov':
        modularity = []
        for thresh in thresholded_arrays:
            print(thresh)
            G, _ = networx(thresh, nodeLabel=None)
            _, mc_clusters = algorithms.markov(G) #use the mc library that is in algorithms
            modularity.append(nx.algorithms.community.modularity(G, mc_clusters))
    else:
        modularity = []
    return percentiles, modularity


def combine_node_attrs(node_attrs_df, WMDz_PC_df, Allens, glob_eff):
    final_df = pd.merge(node_attrs_df, WMDz_PC_df, left_index=True,
                        right_index=True)  # You need the two dfs from Results & Hubs functions
    final_df['Allen_ROI'] = Allens  # This is the list that comes from the unpickled Allen_ROI_dict
    final_df['Delta_Global_Efficiency'] = glob_eff
    # reorder all of the columns to your liking
    final_df = final_df[
        ["Allen_ROI", "Degree", "Betweenness", "Eigenvector_Centrality", "Closeness", "Clustering_Coefficient",
         "WMDz", "PC", "Delta_Global_Efficiency", "Hub_Score"]]
    return final_df


def node_attrs_to_csv(final_df, folder, var_name):
    final_df.to_csv(folder + '/' + var_name + '.csv')
    return


'''def gephiMyNetwork(threshold_matrix, path):
    # I don't know why but gephi function only works when we make the network through here, but everything is the same.
    Gg = nx.to_networkx_graph(threshold_matrix)
    nx.write_gexf(Gg, path)'''