"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

import random

import markov_clustering as mc
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.special as sc
import seaborn as sns
import statsmodels.api as sm
# import necessary libraries
from numpy.random import random
from sklearn.cluster import AgglomerativeClustering
from statsmodels.graphics.regressionplots import abline_plot
from statsmodels.sandbox.stats.multicomp import multipletests


# simple function for loading our csv file
def loadData(data):
    data = pd.read_csv(data)
    data = data.apply(lambda x: x.fillna(x.mean()), axis=0)
    node_names = data.columns.to_list()
    node_number = list(item for item in range(0, len(node_names)))
    nodes = {node_number[i]: node_names[i] for i in range(len(node_number))}
    return data, nodes


# correlate our c-Fos counts between brain regions, df for data
# type for correlation coefficient i.e. "pearson"
def corrMatrix(data):
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
    return rVal, p, p_adjusted, alpha_corrected


# using this function we will threshold based off of p-values previously calculated
def significanceCheck(p_adjusted, corr, alpha, threshold=0.0, names=None, plot=False, include_Negs=True):
    p_adjusted = np.where((p_adjusted >= alpha), 0, p_adjusted)  # if not significant --> zero
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
        else:
            pandas_matrix = pd.DataFrame(threshold_matrix)
            sns.clustermap(pandas_matrix)
            return threshold_matrix, pandas_matrix
    else:
        return threshold_matrix


# pandas_matrix = pd.DataFrame(threshold_matrix, index=list(names.values()), columns=list(names.values()))
# sns.color_palette("viridis", as_cmap=True) sns.clustermap(pandas_matrix,cmap = "viridis",method='ward',
# metric='euclidean',figsize=(10,10),cbar_pos=(.9,.9,.02,.10)) fig1 = plt.figure() sns.heatmap(pandas_matrix,
# cmap = "viridis",square=True,xticklabels=True,yticklabels=True)


def hierarch_clust(threshold_matrix, nodes, allen_groups, plot=False):
    # Do some cute for-loop to examine a variety of distances to cut along the dendrogram and examine the changes in
    # cluster # across cuts
    distances = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    num_clusts = []
    mods = []
    for i in distances:
        hc = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=i, compute_distances=True,
                                     compute_full_tree=True)
        org_hc = hc.fit_predict(threshold_matrix)
        mapped_nodes = {node: cluster for node, cluster in zip(nodes.values(), org_hc)}
        sorted_mapped_nodes = {node: cluster for node, cluster in
                               sorted(mapped_nodes.items(), key=lambda item: item[1])}
        list_of_node_sets = [{node for node, cluster in sorted_mapped_nodes.items() if cluster == j} for j in
                             range(0, max(sorted_mapped_nodes.values()) + 1)]
        mods.append(nx.algorithms.community.modularity(G, list_of_node_sets))
        print(org_hc)  # Here is actually where you make the rules for the clustering
        print(np.bincount(org_hc))  # output to indicate how many nodes fall in each clusters
        num_clusters = len(np.unique(org_hc))
        num_clusts.append(num_clusters)
    modularities = {dist: mod for dist, mod in zip(distances, mods)}
    max_mod_dist = max(modularities, key=modularities.get)
    d = {"Distance": distances, "Number of Clusters": num_clusts}
    df_clusts = pd.DataFrame(d, columns=["Distance", "Number of Clusters"])
    # Next, perform the HC on the distance that is "halfway" down the dendrogram
    hc_2 = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=max_mod_dist, compute_distances=True,
                                   compute_full_tree=True)
    org_hc_2 = hc_2.fit_predict(threshold_matrix)
    nodes_items = nodes.items()  # Now we conduct some tomfoolery to identify clusters in nodes
    nodes_list = list(nodes_items)
    nodes_df = pd.DataFrame(nodes_list)
    nodes_df["cluster"] = org_hc_2
    nodes_df["Allen Group Name"] = allen_groups
    if plot:
        fig2 = plt.figure()
        sch.dendrogram(sch.linkage(threshold_matrix, method='ward'))
        fig3 = plt.figure()
    return df_clusts, nodes_df


# we will create our undirected network graphs based on our matrices
def networx(corr_data, nodeLabel):
    # corr_data = np.arctanh(corr_data) # Another Fisher transformation
    # will zero out the weaker edge connections and also not look at negative edge connections
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)  # use the updated corr_data to make a graph
    graph = nx.relabel_nodes(graph, nodeLabel)
    remove = [node for node, degree in graph.degree() if degree < 1]
    graph.remove_nodes_from(remove)
    pos = nx.spring_layout(graph)
    nx.set_node_attributes(graph, pos, name='pos')
    return graph, pos


def markov(graph, plot=False):
    matrix = nx.to_scipy_sparse_matrix(graph)  # Will generate an adjacency matrix from the graph
    inflation_values = []
    modularity_values = []
    for inflation in [i / 10 for i in range(15, 135, 5)]:
        result = mc.run_mcl(matrix, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        inflation_values.append(inflation)
        modularity_values.append(Q)
    d = {"Inflation": inflation_values, "Modularity": modularity_values}
    df = pd.DataFrame(d, columns=["Inflation", "Modularity"])  # Make a df of the inflation and modularity values
    column = df["Modularity"]
    max_index = column.idxmax()
    optimal_inflation = df["Inflation"].iloc[max_index]
    mc_results = mc.run_mcl(matrix, inflation=optimal_inflation)
    mc_clusters = mc.get_clusters(mc_results)
    if plot:
        numnodes = graph.number_of_nodes()
        positions = {i: (random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
        mc.draw_graph(matrix, mc_clusters, pos=positions, node_size=100, with_labels=True, edge_color='silver')
    return df, mc_clusters


def grab_color_attributes(cluster_list, node_dict):
    color_list = [color for color in sns.color_palette('colorblind', len(cluster_list))]
    color_dict = {}
    for i in range(0, len(cluster_list)):
        for j in cluster_list[i]:
            color_dict[node_dict[j]] = color_list[i]
    color_dict_sorted = {area: color for area, color in sorted(color_dict.items(), key=lambda ele: ele[0])}
    return color_dict_sorted


def grab_attributes(graph):
    deg = nx.degree_centrality(graph)
    between = nx.betweenness_centrality(graph)
    eig = nx.eigenvector_centrality(graph)
    deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[1])}
    between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[1])}
    eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[1])}
    return deg_sort, between_sort, eig_sort


def graph_network(G, color_list):
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    pos = nx.spring_layout(G, k=0.5, seed=3847897236)
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=edge_colors)
    nx.draw_networkx_nodes(G, pos=pos, node_size=node_sizes, node_color=color_list)
    nx.draw_networkx_labels(G, pos=pos)
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.show()
    plt.axis('off')
    return


def get_ordered_degree_list(G):
    degree_ordered = {k: v for k, v in sorted(dict(G.degree()).items(), key=lambda item: item[1])}
    return degree_ordered


def in_silico_deletion(G, plot=False):
    degree_list = [degree for degree in dict(G.degree).values()]
    og_global_eff = nx.global_efficiency(G)
    delta_global_eff = [abs(nx.global_efficiency(nx.from_pandas_adjacency(disruptPropagate(G, node))) - og_global_eff)
                        for node in list(G.nodes())]
    degree_list_const = sm.tools.add_constant(degree_list)
    my_model = sm.OLS(delta_global_eff, degree_list_const).fit()
    print(my_model.summary())
    if plot:
        sns.set()
        fig, ax = plt.subplots()
        plt.scatter(degree_list, delta_global_eff)
        abline_plot(model_results=my_model, ax=ax)
    return delta_global_eff


# this is the disruption propagation model from Vetere et al. 2018
def disruptPropagate(G, target):
    G = nx.to_pandas_adjacency(G)  # create df with node names
    # keep track of negative weights, create a "negs" matrix that encodes where negative  correlations were in the
    # correlation matrix
    Negs = G
    Negs = np.where((Negs < 0), -1, Negs)
    Negs = np.where((Negs > 0), 1, Negs)
    # keep track of edge updates
    Update = np.ones(G.shape)
    np.fill_diagonal(Update, 0)
    # make a copy of the matrix
    G_0 = abs(G.copy())
    # remove edges on target node
    G.loc[target, :] = 0
    G.loc[:, target] = 0
    deltaEdge = pd.DataFrame(G_0 - G)  # edge change df
    Update = np.where((deltaEdge > 0), 0, Update)
    iterG = [G_0, G]  # list that stores iterations of disrupt propagate model
    while True:
        xGraph = iterG[-1]  # arbitrary variable that stores the latest iteration of the model
        yGraph = deltaEdge.tail(len(deltaEdge))  # arbitrary variable that stores the latest entry in the deltaEdge df
        # calculate edge weights
        deltaEdge_sum = yGraph.sum()
        deltaEdge_sum = deltaEdge_sum.to_numpy()
        edgeWeights = np.add.outer(deltaEdge_sum, deltaEdge_sum)
        # calculate node weights
        nodeWeights = xGraph.sum()
        nodeWeights = nodeWeights.to_numpy()
        nodeWeights = np.add.outer(nodeWeights, nodeWeights)
        nodeWeights = nodeWeights - (2 * xGraph)
        # update edge weights
        UpdateMult = Update
        UpdateMult = np.where((Update > 1), 1, UpdateMult)
        UpdateMult = np.where((Update < 0), 0, UpdateMult)

        deltaEdge_x = edgeWeights / nodeWeights
        deltaEdge_x.replace(np.nan, 0)
        deltaEdge_x.replace([np.inf, -np.inf], 0, inplace=True)
        deltaEdge_x = deltaEdge_x * xGraph * UpdateMult

        xGraph = xGraph - deltaEdge_x
        xGraph = np.where((xGraph < 0), 0, xGraph)
        # update matrices
        UpdateSub = Update * 0
        UpdateSub = np.where((deltaEdge_x > 0), 1, UpdateSub)
        Update = Update - UpdateSub
        Update = np.where((xGraph == 0), 0, Update)
        iterG.append(pd.DataFrame(xGraph))
        deltaEdge.append(yGraph)
        if iterG[-1].equals(iterG[-2]):
            break
        finalMat = np.multiply(Negs, iterG[-1])
        finalMat = pd.DataFrame(finalMat)

    return finalMat


# THIS IS FOR TEST PURPOSES ONLY
file = '/home/ryansenne/PycharmProjects/Networks/ChR2_Large_Network.csv'
file2 = '/home/ryansenne/PycharmProjects/Networks/Control_Large_Network.csv'
allen_groups = pd.read_csv('/home/ryansenne/PycharmProjects/Networks/ROIs.csv')

test_data, test_nodes = loadData(file2)
rvals, p, p_adj, rej = corrMatrix(test_data)
threshold_matrix = significanceCheck(p_adj, rvals, 0.001, names=test_nodes, include_Negs=False)
G, pos = networx(threshold_matrix, test_nodes)
# my_del = in_silico_deletion(G, plot=True)
my_list = get_ordered_degree_list(G)
clust, result = hierarch_clust(threshold_matrix, test_nodes, allen_groups['Allen Group Name'])
