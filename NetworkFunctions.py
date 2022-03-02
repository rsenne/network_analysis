"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

# import necessary libraries
from numpy.random import random
import networkx as nx
import pandas as pd
import numpy as np
import scipy.special as sc
import seaborn as sns
import matplotlib.pyplot as plt
import markov_clustering as mc
from bokeh.io import output_file, show
from bokeh.models import (Circle, MultiLine, Plot, Range1d, BoxZoomTool, ResetTool, PanTool)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx
from bokeh.models import ColumnDataSource, LabelSet
from IPython import get_ipython

# THIS IS FOR TEST PURPOSES ONLY
file = '/home/ryansenne/PycharmProjects/Networks/ChR2_Large_Network.csv'
file2 = '/home/ryansenne/PycharmProjects/Networks/Control_Large_Network.csv'
get_ipython().run_line_magic('matplotlib', 'qt')


# simple function for loading our csv file
def loadData(data):
    data = pd.read_csv(data)
    # data.drop(columns=data.columns[data.nunique() <= 3], inplace=True)
    data = data.apply(lambda x: x.fillna(x.mean()), axis=0)
    node_names = data.columns.to_list()
    node_number = list(item for item in range(0, len(node_names)))
    nodes = {node_number[i]: node_names[i] for i in range(len(node_number))}
    return data, nodes


# correlate our c-Fos counts between brain regions, df for data
# type for correlation coefficient i.e. "pearson"
def corrMatrix(data):
    rVal = np.corrcoef(data, rowvar=False)  # calculate pearson coefficients
    rf = rVal[np.triu_indices(rVal.shape[0], 1)]  # upper triangular matrix of data to shuffle for p-value calc
    df = data.shape[1] - 2  # calculate degrees of freedom
    ts = rf * rf * (df / (1 - rf * rf))  # calculate t's
    pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))  # calculate p's from beta incomplete function
    # generate p-value matrix
    p = np.zeros(shape=rVal.shape)
    p[np.triu_indices(p.shape[0], 1)] = pf
    p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]
    p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
    return rVal, p


# using this function we will threshold based off of p-values previously calculated
def significanceCheck(p, corr, alpha, threshold=0.0, names=None, plot=False, include_Negs=True):
    p = np.where((p >= alpha), 0, p)  # if not significant --> zero
    p = np.where((p != 0), 1, p)  # if significant --> one
    if not include_Negs:
        p = np.where((corr < 0), 0, p)
    threshold_matrix = np.multiply(corr, p)  # remove any insignificant correlations
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


# we will create our undirected network graphs based on our matrices
def networx(corr_data, nodeLabel):
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)
    graph = nx.relabel_nodes(graph, nodeLabel)
    pos = nx.random_layout(graph)
    nx.set_node_attributes(graph, pos, name='pos')
    return graph, pos


def grab_attributes(graph):
    deg = nx.degree_centrality(graph)
    between = nx.betweenness_centrality(graph)
    eig = nx.eigenvector_centrality(graph)
    deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[1])}
    between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[1])}
    eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[1])}
    return deg_sort, between_sort, eig_sort

def graph_network(G):
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    pos = nx.spring_layout(G, k=0.5, seed=3847897236)
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos, width=1, edge_color=edge_colors)
    nx.draw_networkx_nodes(G, pos=pos, node_size=node_sizes)
    nx.draw_networkx_labels(G, pos=pos)
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.show()
    plt.axis('off')
    return



# # calculate the necessary graph attributes such as centrality, betweenness, global efficiency etc.
# def graphAttributes(graph, target, threshold, iterations, mode):


# this is the disruption propagation model from Vetere et al. 2018
def disruptPropagate(G, target, nodeNames):
    G = pd.DataFrame(G, index=list(nodeNames.values()), columns=list(nodeNames.values()))  # create df with node names
    # keep track of negative weights, create a "negs" matrix that encodes where negative correlations were in the
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


# this is all for testing please ignore it
# a, b = loadData(file)
# c, d = corrMatrix(a)
# e, m = significanceCheck(d, c, 0.0001, 0.0, b, True, True)
# q, r = networx(e, b)
# yo = GraphingNetwork(q, "Chr2", nx.circular_layout)

aa, bb = loadData(file2)
cc, dd = corrMatrix(aa)
ee, mm = significanceCheck(dd, cc, 0.0001, 0.5, bb, True, True)
qq, rr = networx(ee, bb)
yoyo = graph_network(qq)

# dxm = e - ee
# numnodes=qq.number_of_nodes()
# positions = {i:(random() * 2 - 1, random() * 2 - 1) for i in range(numnodes)}
# matrix = nx.to_scipy_sparse_matrix(qq)
# result = mc.run_mcl(matrix, inflation=6)           # run MCL with default parameters
# clusters = mc.get_clusters(result)
# Q = mc.modularity(matrix=result, clusters=clusters)
# print("inflation:", 20, "modularity:", Q)
# print("inflation:", inflation, "modularity:", Q)
# plt.figure()
# mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")

# inflation_values = []
# for inflation in [i for i in np.linspace(6, 7, 10)]:
#     result = mc.run_mcl(matrix, inflation=inflation)
#     clusters = mc.get_clusters(result)
#     Q = mc.modularity(matrix=result, clusters=clusters)
#     inflation_values.append(Q)
#     print("inflation:", inflation, "modularity:", Q)

# plt.plot(Q)


#   dd = [[1, 2, 3],
#         [4, 5, 6],
#         [7, 8, 9]]
# GraphingNetwork(dd)
# h = disruptPropagate(dd, 0)

# tyty = disruptPropagate(e, "PAG", b)
# mn, nm = networx(tyty.to_numpy(), b)
# # graphingNetwork(mn)
