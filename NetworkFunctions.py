"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

# import necessary libraries
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
from numpy.random import random
from sklearn.cluster import AgglomerativeClustering
from statsmodels.graphics.regressionplots import abline_plot
from statsmodels.sandbox.stats.multicomp import multipletests
import pickle as pkl
import matplotlib.patches as mpatches

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
def significanceCheck(p_adjusted, corr, alpha, threshold=0.0, names=None, plot=False, include_Negs=True, Anatomy=None):
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
            if Anatomy:
                #Create sorted list from the unpickled dictionary
                sorted_dict = dict(sorted(Anatomy.items(),key=lambda item: item[1]))
                list_for_sort = list(sorted_dict.keys())

                allen_pandas = pandas_matrix[list_for_sort].loc[list_for_sort] #new pandas matrix that is organized by allen

                #plot a different correlation matrix
                allens = list(sorted_dict.values())
                allens_unique = np.unique(allens)
                color_list = [color for color in sns.color_palette('Set3',n_colors = len(allens_unique))]
                color_ref = dict(zip(map(str, allens_unique),color_list))
                allen_colors = pd.Series(allens,index = allen_pandas.columns).map(color_ref)

                #create legends for Allen groups with the colors from the color_ref dictionary
                cerebellum = mpatches.Patch(color=color_list[0],label='Cerebellum')
                cort_plate = mpatches.Patch(color=color_list[1],label='Cortical Plate')
                cort_subplate = mpatches.Patch(color=color_list[2],label='Cortical Subplate')
                hypothalamus = mpatches.Patch(color=color_list[3],label='Hypothalamus')
                medulla = mpatches.Patch(color=color_list[4],label='Medulla')
                midbrain = mpatches.Patch(color=color_list[5],label='Midbrain')
                pallidum = mpatches.Patch(color=color_list[6],label='Pallidum')
                pons = mpatches.Patch(color=color_list[7],label='Pons')
                striatum = mpatches.Patch(color=color_list[8],label='Striatum')
                thal = mpatches.Patch(color=color_list[9],label='Thalamus')

                # plot the new corrleation matrix
                plt.figure(1)
                sns.clustermap(allen_pandas, cmap='viridis', row_colors=allen_colors, col_colors=allen_colors,
                               row_cluster=False, col_cluster=False, xticklabels=False, yticklabels=False,
                               figsize=(10, 10),cbar_pos=(0.1,0.15,.02,.4),cbar_kws={'label':'Pearson Correlation (R)'})
                plt.legend(handles=[cerebellum, cort_plate, cort_subplate, hypothalamus, medulla,
                                    midbrain, pallidum, pons, striatum, thal],bbox_to_anchor =(5.0,1.6))
        else:
            pandas_matrix = pd.DataFrame(threshold_matrix)
        fig = sns.clustermap(pandas_matrix, cmap='viridis', method='ward', metric='euclidean', figsize=(10,10), cbar_pos=(.9,.9,.02,.10))
        return threshold_matrix, pandas_matrix
    else:
        return threshold_matrix


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


def hierarch_clust(G, nodes, allen_groups, plot=False):
    #Create the adjacency matrix from the imported graph
    adj_array = nx.to_numpy_array(G)
    # cluster # across cuts
    distances = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65] #list of distances to cut along the dendrogram
    num_clusts = []
    mods = []
    #create a for loop to iterate through the distances with HC and append to the two empty lists above
    for i in distances:
        hc = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=i, compute_distances=True,
                                     compute_full_tree=True)
        org_hc = hc.fit_predict(adj_array)
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
    hc_2 = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=max_mod_dist,
                                   compute_distances=True,
                                   compute_full_tree=True)
    org_hc_2 = hc_2.fit_predict(adj_array)
    nodes_items = nodes.items()  # Now we conduct some tomfoolery to identify clusters in nodes
    nodes_list = list(nodes_items)
    nodes_df = pd.DataFrame(nodes_list)
    nodes_df["cluster"] = org_hc_2
    nodes_df["Allen Group Name"] = allen_groups
    if plot:
        plt.figure(1)
        sch.dendrogram(sch.linkage(adj_array, method='ward'))
    return df_clusts, nodes_df


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

def louvain(graph,plot=False):
    


def grab_color_attributes(cluster_list, node_dict):
    color_list = [color for color in sns.color_palette('colorblind', len(cluster_list))]
    color_dict = {}
    for i in range(0, len(cluster_list)):
        for j in cluster_list[i]:
            color_dict[node_dict[j]] = color_list[i]
    color_dict_sorted = {area: color for area, color in sorted(color_dict.items(), key=lambda ele: ele[0])}
    return color_dict_sorted


def grab_attributes(graph):
    #attributes for nodes
    deg = nx.degree_centrality(graph)
    between = nx.betweenness_centrality(graph)
    eig = nx.eigenvector_centrality(graph)
    deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[1])}
    between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[1])}
    eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[1])}
    #attributes for the whole graph
    avrg_clust_coeff = nx.average_clustering(graph,weight='weight')
    global_eff = nx.global_efficiency(graph)
    return deg_sort, between_sort, eig_sort,avrg_clust_coeff,global_eff


def graph_network(G, color_list, pos_dict):
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    # pos = nx.spring_layout(G, k=0.5, seed=3847897236)
    # print(pos)
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos_dict, width=1, edge_color=edge_colors, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(G, pos=pos_dict, node_size=node_sizes, node_color=color_list, linewidths=1,
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos=pos_dict)
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

def get_position_data(cluster_list, node_names):
    number_of_clusters = len(cluster_list)
    nodes_list = [x for x in range(0, number_of_clusters)]
    pos_graph = nx.Graph()
    pos_graph.add_nodes_from(nodes_list)
    pos = nx.circular_layout(pos_graph, scale=30, dim=2)
    num_of_nodes = [len(node) for node in cluster_list]
    point_clouds = [get_point_cloud(lens) for lens in num_of_nodes]
    for i in range(len(point_clouds)):
        for j in range(len(point_clouds[i])):
            point_clouds[i][j][0] += pos[i][0]
            point_clouds[i][j][1] += pos[i][1]
    point_cloud_map = {cluster: pos_list for cluster, pos_list in enumerate(point_clouds)}
    pos_dict = {}
    for i in range(len(cluster_list)):
        for j in range(len(cluster_list[i])):
            pos_dict.update({node_names[cluster_list[i][j]]: np.array(point_cloud_map[i][j])})
    nx.rescale_layout_dict(pos_dict, 10)
    return pos_dict


def sunflower_theta(n):
    golden_ratio = ((1 + 5 ** 0.5) / 2) ** 2
    return 2 * np.pi / golden_ratio * n


def sunflower_r(n, c=1):
    return c * (n ** 0.5)


def get_point_cloud(k=0):
    n = [i for i in range(1, k + 1)]
    r = [sunflower_r(i) for i in n]
    theta = [sunflower_theta(j) for j in n]
    point_cloud_x = (r * np.cos(theta))
    point_cloud_y = (r * np.sin(theta))
    point_cloud = [[x, y] for x, y in zip(point_cloud_x, point_cloud_y)]
    return point_cloud


# THIS IS FOR TEST PURPOSES ONLY

file = '/home/ryansenne/PycharmProjects/Networks/ChR2_Large_Network.csv'
file2 = '/home/ryansenne/PycharmProjects/Networks/Control_Large_Network.csv'
allen_groups = pd.read_csv('/home/ryansenne/PycharmProjects/Networks/ROIs.csv')

test_data, test_nodes = loadData(file2)
rvals, p, p_adj, rej = corrMatrix(test_data)
threshold_matrix = significanceCheck(p_adj, rvals, 0.001, names=test_nodes, include_Negs=False)
G, pos = networx(threshold_matrix, test_nodes)
df, mark_clust = markov(G)
color_list = grab_color_attributes(mark_clust, test_nodes)
pos_dict = get_position_data(mark_clust, test_nodes)
graph_network(G, list(color_list.values()), pos_dict)
# my_del = in_silico_deletion(G, plot=True)
my_list = get_ordered_degree_list(G)
clust, result = hierarch_clust(threshold_matrix, test_nodes, allen_groups['Allen Group Name'])'''
# clust, result = hierarch_clust(threshold_matrix, test_nodes, allen_groups['Allen Group Name'])
