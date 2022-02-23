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
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import markov_clustering as mc
from bokeh.io import output_file, show
from bokeh.models import (Circle, MultiLine, Plot, Range1d, BoxZoomTool, ResetTool)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx
from bokeh.models import ColumnDataSource, LabelSet
import random
from IPython import get_ipython

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
    rVal = np.corrcoef(data, rowvar=False) # calculate pearson coefficients
    rVal[np.isnan(rVal)] = 0 # Will make all NaN values into zero
    rf = rVal[np.triu_indices(rVal.shape[0], 1)]  # upper triangular matrix of data to shuffle for p-value calc
    df = data.shape[1] - 2  # calculate degrees of freedom
    ts = rf * rf * (df / (1 - rf * rf))  # calculate t's
    pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))  # calculate p's from beta incomplete function
    # generate p-value matrix
    p = np.zeros(shape=rVal.shape)
    p[np.triu_indices(p.shape[0], 1)] = pf
    p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]
    p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
    #Multiple comparison of p values using Bonferroni correction
    rejected, p_adjusted, _, alpha_corrected = multipletests(p, alpha = 0.05, method= 'bonferroni', is_sorted = True)
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
        pandas_matrix = pd.DataFrame(threshold_matrix, index=list(names.values()), columns=list(names.values()))
        sns.color_palette("viridis", as_cmap=True)
        sns.clustermap(pandas_matrix,cmap = "viridis",method='ward',metric='euclidean',figsize=(10,10),cbar_pos=(.9,.9,.02,.10))
    return threshold_matrix, pandas_matrix


#I think we can perform hierarchical clustering here?
def hierarch_clust(threshold_matrix,nodes,plot = False):
    #threshold_matrix = np.arctanh(threshold_matrix) # Will do a Fisher z transformation of the data
    #Do some cute for-loop to examine a variety of distances to cut along the dendrogram and examine the changes in cluster # across cuts
    distances = [10,20,30,40,50,60]
    num_clusts = []
    for i in distances:
        hc = AgglomerativeClustering(n_clusters=None,linkage='ward',distance_threshold=i,compute_distances=True,compute_full_tree=True)
        org_hc = hc.fit_predict(threshold_matrix) # Here is actually where you make the rules for the clustering
        print(np.bincount(org_hc)) #output to indicate how many nodes fall in each clustes
        num_clusters = len(np.unique(org_hc))
        num_clusts.append(num_clusters)
    d = {"Distance":distances,"Number of Clusters":num_clusts}
    df_clusts = pd.DataFrame(d,columns=["Distance","Number of Clusters"])
    #Next, perform the HC on the distance that is "halfway" down the dendrogram
    hc_2 = AgglomerativeClustering(n_clusters=None,linkage='ward',distance_threshold=40,compute_distances=True,compute_full_tree=True)
    org_hc_2 = hc_2.fit_predict(threshold_matrix)
    nodes_items = nodes.items() #Now we conduct some tomfoolery to identify clusters in nodes
    nodes_list = list(nodes_items)
    nodes_df = pd.DataFrame(nodes_list)
    nodes_df["cluster"] = org_hc_2
    #Reduce the data to two dimensions using PCA
    pca = PCA(n_components=2)
    components = pca.fit_transform(threshold_matrix)
    if plot:
        fig1 = plt.figure()
        sch.dendrogram(sch.linkage(threshold_matrix, method='ward'))
        fig2 = plt.figure()
        plt.scatter(x = components[:,0], y = components[:,1],c=nodes_df["cluster"],cmap = "rainbow")
    return df_clusts, nodes_df, components


# we will create our undirected network graphs based on our matrices
def networx(corr_data,nodeLabel,r_thresh = 0.85):
    #corr_data = np.arctanh(corr_data) # Another Fisher transformation
    corr_data = np.where((corr_data < r_thresh),0,corr_data) #will zero out the weaker edge connections and also not look at negative edge connections
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph) #use the updated corr_data to make a graph
    graph = nx.relabel_nodes(graph, nodeLabel)
    pos = nx.spring_layout(graph)
    nx.set_node_attributes(graph, pos, name='pos')
    return graph, pos


def markov(graph,plot = False):
    numnodes = graph.number_of_nodes()
    positions = {i: (random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
    matrix = nx.to_scipy_sparse_matrix(graph) #Will generate an adjacency matrix from the graph
    inflation_values = []
    modularity_values = []
    for inflation in [i /10 for i in range(15,135,5)]:
        result = mc.run_mcl(matrix, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        inflation_values.append(inflation)
        modularity_values.append(Q)
    d = {"Inflation":inflation_values,"Modularity":modularity_values}
    df = pd.DataFrame(d,columns=["Inflation","Modularity"]) #Make a df of the inflation and modularity values
    column = df["Modularity"]
    max_index = column.idxmax()
    optimal_inflation = df["Inflation"].iloc[max_index]
    mc_results = mc.run_mcl(matrix,inflation = optimal_inflation)
    mc_clusters = mc.get_clusters(mc_results)
    if plot:
        mc.draw_graph(matrix,mc_clusters,pos = positions,node_size = 100,with_labels = True,edge_color = 'silver')
    return df,mc_results,mc_clusters


def GraphingNetwork(G, plot_title, nx_layout):
    negativeCorr, positiveCorr = 'red', 'black'
    edge_attribs = {}
    for i, j, _ in G.edges(data=True):
        edge_color = negativeCorr if G[i][j]['weight'] < 0 else positiveCorr
        edge_attribs[(i, j)] = edge_color
    nx.set_edge_attributes(G, edge_attribs, "edge_color")
    plot = Plot(plot_width=2000, plot_height=2000,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))
    plot.title.text = plot_title
    plot.add_tools(BoxZoomTool(), ResetTool())
    graph_renderer = from_networkx(G, nx_layout, scale=1, center=(0, 0))
    x, y = zip(*graph_renderer.layout_provider.graph_layout.values())
    source = ColumnDataSource({'x': pd.Series(x), 'y': pd.Series(y), 'field': list(G.nodes())})
    labels = LabelSet(x='x', y='y', text_font_size='14px', x_offset=0, y_offset=0, text='field', source=source)
    graph_renderer.node_renderer.glyph = Circle(size=50, fill_color=Spectral4[0])
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=1)
    plot.renderers.append(graph_renderer)
    plot.renderers.append(labels)
    output_file("interactive_graphs.html")
    show(plot)

# # calculate the necessary graph attributes such as centrality, betweenness, global efficiency etc.
# def graphAttributes(graph, target, threshold, iterations, mode):
# For any sort of quantitative statistics, need to Fisher transform the correlation data


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

# THIS IS FOR TEST PURPOSES ONLY
#file = 'TempChR2.csv'
#file2 = 'TempControl.csv'
#get_ipython().run_line_magic('matplotlib', 'qt')

# this is all for testing please ignore it
# a, b = loadData(file)
# c, d = corrMatrix(a)
# e, m = significanceCheck(d, c, 0.0001, 0.0, b, True, True)
# q, r = networx(e, b)
# yo = GraphingNetwork(q, "Chr2", nx.circular_layout)

#aa, bb = loadData(file2)
#cc, dd = corrMatrix(aa)
#ee, mm = significanceCheck(dd, cc, 0.0001, 0.0, bb, True, True)
#qq, rr = networx(ee, bb)
# yoyo = GraphingNetwork(qq, "control", nx.circular_layout)

# dxm = e - ee
#numnodes=qq.number_of_nodes()
#ChR2_
#matrix = nx.to_scipy_sparse_matrix(qq)
#result = mc.run_mcl(matrix, inflation=6)           # run MCL with default parameters
#clusters = mc.get_clusters(result)
# Q = mc.modularity(matrix=result, clusters=clusters)
# print("inflation:", 20, "modularity:", Q)
#print("inflation:", inflation, "modularity:", Q)
# plt.figure()
# mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")



#plt.plot(Q)


#   dd = [[1, 2, 3],
#         [4, 5, 6],
#         [7, 8, 9]]
# GraphingNetwork(dd)
# h = disruptPropagate(dd, 0)

# tyty = disruptPropagate(e, "PAG", b)
# mn, nm = networx(tyty.to_numpy(), b)
# # graphingNetwork(mn)
