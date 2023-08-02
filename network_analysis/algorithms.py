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
<<<<<<< Updated upstream
import os

=======
from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse import csr_matrix
from numpy import square
from numpy import trace
from numpy import amax
from math import sqrt
>>>>>>> Stashed changes



def hierarch_clust(graph, nodes, allen_groups, plot=False):
    adj_matrix = nx.to_numpy_matrix(graph)
    distances = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
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
    if plot:
        plt.figure()
        sch.dendrogram(sch.linkage(adj_matrix,method='ward',metric='euclidean'))
    return df_clust_cuts,clust_assigns,sorted(list(clusters))


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
<<<<<<< Updated upstream
    if plot:
        numnodes = graph.number_of_nodes()
        positions = {i: (random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
        mc.draw_graph(matrix, mc_clusters, pos=positions, node_size=100, with_labels=True, edge_color='silver')
    return df, mc_clusters
=======

    # Generate an N x 1 vector of cluster ids
    node_nums = list(nodes.keys())
    mc_clust_index = [i for i in range(len(mc_clusters))]
    cluster_vector_mc = []
    for node in node_nums:
        for clust in mc_clust_index:
            if node in mc_clusters[clust]:
                cluster_vector_mc.append(clust)
            else:
                pass
    clust_dict = {list(graph.nodes())[i]: {'cluster': cluster_vector_mc[i]} for i in range(len(cluster_vector_mc))}
    nx.set_node_attributes(graph, clust_dict)
    cluster_vector_mc = np.array(cluster_vector_mc)
    cluster_vector_mc = np.reshape(cluster_vector_mc, (len(cluster_vector_mc), 1))
    # plt.plot(inflation_values, modularity_values)
    return df, mc_clusters, cluster_vector_mc
>>>>>>> Stashed changes


def louvain(graph,nodes):
    node_nums = {value:key for key,value in nodes.items()}
    graph = nx.relabel_nodes(graph,node_nums)
    resolutions = [0.5,1.0,1.2,1.4,1.6,1.8,2.0]
    lou_mod = []
    for i in resolutions:
        lou_clust = nx_comm.louvain_communities(graph,resolution=i)
        lou_mod.append(nx_comm.modularity(graph,lou_clust))
    lou_modularities = {res:mod for res, mod in zip(resolutions,lou_mod)}
    max_res = max(lou_modularities,key=lou_modularities.get)
    max_mod_lou_comm = nx_comm.louvain_communities(graph,resolution=max_res,randomize=True, seed=100)
    max_mod_lou_comm = [tuple(c) for c in max_mod_lou_comm]
    return max_mod_lou_comm


def in_silico_deletion(G, title = '', results_folder = '', plot=False, save = False):
    degree_list = [degree for degree in dict(G.degree).values()]
    og_global_eff = nx.global_efficiency(G)
    delta_global_eff = [abs(nx.global_efficiency(nx.from_pandas_adjacency(disruptPropagate(G, node))) - og_global_eff)
                        for node in list(G.nodes())]
    degree_list_const = sm.tools.add_constant(degree_list)
    my_model = sm.OLS(delta_global_eff, degree_list_const).fit()
    print(title)
    print(my_model.summary())
    if plot:
        sns.set()
        fig, ax = plt.subplots()
        plt.scatter(degree_list, delta_global_eff)
        abline_plot(model_results=my_model, ax=ax)
        plt.title(title, fontsize=20)
        if save:
            fig.savefig(os.path.join(results_folder, 'Insilico deletion' + title +\
                                     '.png'), dpi=300)
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
<<<<<<< Updated upstream
=======


def Similarity(A1, A2):
    '''
    use deltacon0 to compute similarity
    CITATION: Danai Koutra, Joshua T. Vogelstein, Christos Faloutsos
    DELTACON: A Principled Massive-Graph Similarity Function
    '''
    S1 = InverseMatrix(A1)
    S2 = InverseMatrix(A2)

    d = 0
    for i in range(A1.shape[0]):
        for j in range(A1.shape[0]):
            d += (sqrt(S1[i, j]) - sqrt(S2[i, j])) ** 2
    d = sqrt(d)
    sim = 1 / (1 + d)
    return sim


def InverseMatrix(A):
    '''
    use Fast Belief Propagatioin
    CITATION: Danai Koutra, Tai-You Ke, U. Kang, Duen Horng Chau, Hsing-Kuo
    Kenneth Pao, Christos Faloutsos
    Unifying Guilt-by-Association Approaches
    return [I+a*D-c*A]^-1
    '''
    A = csr_matrix(A)
    I = identity(A.shape[0])  # identity matirx
    D = diags(sum(A).toarray(), [0])  # diagonal degree matrix

    c1 = trace(D.toarray()) + 2
    c2 = trace(square(D).toarray()) - 1
    h_h = sqrt((-c1 + sqrt(c1 * c1 + 4 * c2)) / (8 * c2))

    a = 4 * h_h * h_h / (1 - 4 * h_h * h_h)
    c = 2 * h_h / (1 - 4 * h_h * h_h)

    # M=I-c*A+a*D
    # S=inv(M.toarray())
    '''
    compute the inverse of matrix [I+a*D-c*A]
    use the method propose in Unifying Guilt-by-Association equation 5
    '''
    M = c * A - a * D
    S = I
    mat = M
    power = 1
    while amax(M.toarray()) > 10 ** (-9) and power < 7:
        S = S + mat
        mat = mat * M
        power += 1

    return S


def degree_preserving_randomization(G, niter=25000):
    G0 = G.copy()
    i = 0
    while i < niter:
        r1 = np.random.randint(0, 155)
        r2 = np.random.randint(0, 155)
        edge1 = list(G0.edges())[r1]
        edge2 = list(G0.edges())[r2]
        weight1 = G0.edges()[edge1[0], edge1[1]]
        weight2 = G0.edges()[edge2[0], edge2[1]]
        if edge1[0] != edge2[1] and edge1[1] != edge2[0]:
            if G0.has_edge(edge1[0], edge2[1]) or G0.has_edge(edge1[1], edge2[0]):
                continue
            else:
                G0.remove_edge(edge1[0], edge1[1])
                G0.remove_edge(edge2[0], edge2[1])
                G0.add_edge(edge1[0], edge2[1], weight=weight1)
                G0.add_edge(edge1[1], edge2[0], weight=weight2)
                i += 1
        else:
            continue
    return G0
>>>>>>> Stashed changes
