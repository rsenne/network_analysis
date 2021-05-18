"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

# import necessary libraries
import networkx as nx
import pandas as pd
import numpy as np
import scipy.special as sc
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt

# THIS IS FOR TEST PURPOSES ONLY
file = '/home/ryansenne/The_Ramirez_Group/Cell_Counts.csv'


# simple function for loading our csv file
def loadData(data):
    data = pd.read_csv(data)
    node_names = data.columns.to_list()  # these are the names of our nodes, we need these for later
    node_number = list(item for item in range(0, len(node_names)))  # find node number
    nodes = {node_number[i]: node_names[i] for i in range(len(node_number))}  # ordered list of nodes
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
def significanceCheck(p, corr, alpha, threshold=0.0, plot=False):
    p = np.where((p >= alpha), 0, p)  # if not significant --> zero
    p = np.where((p != 0), 1, p)  # if significant --> one
    threshold_matrix = np.multiply(corr, p)  # remove any insignificant correlations
    # remove correlations below threshold
    threshold_matrix = np.where((abs(threshold_matrix) < threshold), 0, threshold_matrix)
    # create a heatmap of correlations if wanted
    if plot:
        sns.heatmap(threshold_matrix, cmap="coolwarm")
        plt.show()
    return threshold_matrix


# we will create our undirected network graphs based on our matrices
def networx(corr_data, nodeLabel):
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)  # create network from a thresholded matrix
    remove = [node for node, degree in dict(graph.degree()).items() if degree == 0]  # find nodes with zero degree
    graph.remove_nodes_from(remove)  # remove all nodes of zero degree
    graph = nx.relabel_nodes(graph, nodeLabel)  # add the names of our nodes from loadData()
    pos = nx.circular_layout(graph)  # change the shape of our graph for plotting, this is arbitrary
    nx.set_node_attributes(graph, pos, name='pos')  # add position data for plotting
    return graph, pos


# this code is taken near-verbatim from the plotly example, I intend to rework the color schemes, but hey, why rebuild
# the wheel?
def graphingNetwork(G):
    global degree_list
    names = []
    for x in G.nodes():
        names.append(x)
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)
        degree_list = []
    for name in names:
        degree_list.append(G.degree(name))
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers + text',
        text=names,
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=np.array(degree_list) * 15,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right',
            ),
            line_width=2))
    node_adjacencies = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
    node_trace.marker.color = node_adjacencies
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='<br>Mouse Chamber ChR2',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.update_layout(font_size=40)
    return fig.show(renderer='browser'), degree_list


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
a, b = loadData(file)
c, d = corrMatrix(a)
e = significanceCheck(d, c, 0.05, 0, True)
q, r = networx(e, b)
yo = graphingNetwork(q)
# god = disruptPropagate()
# dd = [[1, 2, 3],
#       [4, 5, 6],
#       [7, 8, 9]]
tyty = disruptPropagate(e, "PAG", b)
mn, nm = networx(tyty.to_numpy(), b)
graphingNetwork(mn)

