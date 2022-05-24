import networkx as nx
import networkx.algorithms.community as nx_comm
import pandas as pd
import numpy as np
import markov_clustering as mc
import statsmodels.api as sm
import scipy.cluster.hierarchy as sch
import random
import seaborn as sns
import matplotlib.pyplot as plt
from numpy.random import random
from statsmodels.graphics.regressionplots import abline_plot
from sklearn.cluster import AgglomerativeClustering


def hierarch_clust(graph, nodes, allen_groups, plot=False):
    adj_matrix = nx.to_numpy_array(graph) #Convert your graph back into a numpy array
    distances = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65] #set a list of arbitrary distances to cut along dendrogram
    #Run through a series of clustering using the various distance metrics and calculate the max modularity
    num_clusts = []
    mods = []
    for i in distances:
        hc = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=i, compute_distances=True,
                                     compute_full_tree=True)
        org_hc = hc.fit_predict(adj_matrix)
        mapped_nodes = {node: cluster for node, cluster in zip(nodes.values(), org_hc)}
        sorted_mapped_nodes = {node: cluster for node, cluster in
                               sorted(mapped_nodes.items(), key=lambda item: item[1])}
        list_of_node_sets = [{node for node, cluster in sorted_mapped_nodes.items() if cluster == j} for j in
                             range(0, max(sorted_mapped_nodes.values()) + 1)]
        mods.append(nx.algorithms.community.modularity(graph, list_of_node_sets))
        num_clusters = len(np.unique(org_hc))
        num_clusts.append(num_clusters)
    modularities = {dist: mod for dist, mod in zip(distances, mods)}
    max_mod_dist = max(modularities, key=modularities.get)
    d = {"Distance Cut": distances, "Number of Clusters": num_clusts}
    df_clust_cuts = pd.DataFrame(d, columns=["Distance Cut", "Number of Clusters"])

    # Next, perform the HC on the distance that is "halfway" down the dendrogram
    hc_2 = AgglomerativeClustering(n_clusters=None, linkage='ward', distance_threshold=max_mod_dist,
                                   compute_distances=True,compute_full_tree=True)
    org_hc_2 = hc_2.fit_predict(adj_matrix)
    nodes_keys = np.array(list(nodes.keys()))
    clust_assigns = pd.DataFrame(zip(nodes_keys,org_hc_2),columns=["Node Number","Cluster Number"])
    clust_assigns["Allen Group Name"] = allen_groups
    clusters = set()
    for i in list(np.unique(org_hc_2)):
        cluster = tuple(clust_assigns.loc[clust_assigns["Cluster Number"]==i,"Node Number"].tolist())
        clusters.add(cluster)

    #Generate a N x 1 vector of cluster ids for calculating WMDz and PC using functions from bctpy
    clust_vector_hc = np.array(clust_assigns['Cluster Number'])
    clust_vector_hc = np.reshape(clust_vector_hc,(155,1))

    if plot:
        plt.figure()
        sch.dendrogram(sch.linkage(adj_matrix,method='ward',metric='euclidean'))

    return df_clust_cuts,modularities,clust_assigns,sorted(list(clusters)),clust_vector_hc


def markov(graph,nodes):
    matrix = nx.to_scipy_sparse_array(graph)  # Will generate an adjacency matrix from the graph
    #Run the clustering algorithm using a range of inflation values and calculate the max modularity with inflation
    inflation_values = []
    modularity_values = []
    for inflation in [i / 10 for i in range(15, 135, 5)]:
        print(inflation)
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

    #Run the algorithm once more with the inflation used that generated the max mod
    mc_results = mc.run_mcl(matrix, inflation=optimal_inflation)
    mc_clusters = mc.get_clusters(mc_results)

    #Generate a N x 1 vector of cluster ids
    node_nums = list(nodes.keys())
    mc_clust_index = [i for i in range(len(mc_clusters))]
    cluster_vector_mc = []
    for node in node_nums:
        for clust in mc_clust_index:
            if node in mc_clusters[clust]:
                cluster_vector_mc.append(clust)
            else:
                pass
    cluster_vector_mc = np.array(cluster_vector_mc)
    cluster_vector_mc = np.reshape(cluster_vector_mc,(155,1))

    return df, mc_clusters, cluster_vector_mc


def louvain(graph,nodes,n_iters):
    node_nums = {value:key for key,value in nodes.items()}
    graph = nx.relabel_nodes(graph,node_nums)
    resolutions = [0.5,1.0,1.2,1.4,1.6,1.8,2.0] #Use different iterations instead of resolutions and average the 1000 Q values
    lou_mod = []
    for i in resolutions:
        lou_clust = nx_comm.louvain_communities(graph,resolution=i)
        lou_mod.append(nx_comm.modularity(graph,lou_clust))
    lou_modularities = {res:mod for res, mod in zip(resolutions,lou_mod)}
    max_res = max(lou_modularities,key=lou_modularities.get)
    louvain_iters = [nx_comm.louvain_communities(graph,resolution=max_res,seed='random_state') for i in range(n_iters)]
    '''
    max_mod_lou_comm = [tuple(c) for c in max_mod_lou_comm]'''

    ''''#Generate a N x 1 np.array to use for PC and WMDz calculations
    node_nums = list(nodes.keys())
    lou_clust_index = [i for i in range(len())]
    cluster_vector_lou = []
    
    for node in node_nums:
        for clust in lou_clust_index:
            if node in lou_clusters[clust]:
                cluster_vector_lou.append(clust)
            else:
                pass
    cluster_vector_lou = np.array(cluster_vector_lou)
    cluster_vector_lou = np.reshape(cluster_vector_lou,(155,1))
    '''
    return max_mod_lou_comm


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
        plt.xlabel('Degree of Node Deleted')
        plt.ylabel(r'$\Delta$' + ' ' + 'Global Efficiency')
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
