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
from scipy.sparse import identity
from scipy.sparse import diags
from numpy import square
from numpy import trace
from numpy import amax
from math import sqrt


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


def markov(graph, nodes):
    matrix = nx.to_scipy_sparse_array(graph)  # Will generate an adjacency matrix from the graph
    #Run the clustering algorithm using a range of inflation values and calculate the max modularity with inflation
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
    # Use different iterations instead of resolutions and average the 1000 Q values
    louvain_iters = [nx_comm.louvain_communities(graph, resolution=1.0) for i in range(n_iters)]
    lou_mod = []
    for i in louvain_iters:
        mods = nx.algorithms.community.modularity(graph,i,weight='weight',resolution=1)
        lou_mod.append(mods)

    #Combine all of the results into a dataframe
    d = {'Louvain_Results':louvain_iters,'Modularity':lou_mod}
    df = pd.DataFrame(d,columns=['Louvain_Results','Modularity'])
    column = df['Modularity']
    max_mod = column.idxmax()
    comm_max_mod = df['Louvain_Results'].iloc[max_mod]

    #Get the value of the max Q
    lou_max_mod = column.max()

    #Average all of the Q values together
    lou_mod_mean = df['Modularity'].mean()

    #Get a list that is used for plotting utils to cluster based on community structure
    max_mod_lou_comm = [tuple(c) for c in comm_max_mod]

    # Generate a N x 1 np.array to use for PC and WMDz calculations
    node_keys = list(nodes.keys())
    lou_clust_index = [i for i in range(len(max_mod_lou_comm))]
    cluster_vector_lou = []
    for node in node_keys:
        for clust in lou_clust_index:
            if node in max_mod_lou_comm[clust]:
                cluster_vector_lou.append(clust)
            else:
                pass
    cluster_vector_lou = np.array(cluster_vector_lou)
    cluster_vector_lou = np.reshape(cluster_vector_lou, (155, 1))

    return max_mod_lou_comm,lou_max_mod,lou_mod_mean,cluster_vector_lou


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
