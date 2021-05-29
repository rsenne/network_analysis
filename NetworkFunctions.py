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
from bokeh.io import output_file, show
from bokeh.models import (BoxZoomTool, Circle, HoverTool,
                          MultiLine, Plot, Range1d, ResetTool,)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx
from bokeh.models import ColumnDataSource, LabelSet

# THIS IS FOR TEST PURPOSES ONLY
file = 'test.csv'


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
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)
    graph = nx.relabel_nodes(graph, nodeLabel)
    pos = nx.circular_layout(graph)
    nx.set_node_attributes(graph, pos, name='pos')
    return graph, pos


def GraphingNetwork(G):
    negativeCorr, positiveCorr = 'red', 'black'
    edge_attribs = {}
    for i, j, _ in G.edges(data=True):
        edge_color = negativeCorr if G[i][j]['weight'] < 0 else positiveCorr
        edge_attribs[(i, j)] = edge_color
    nx.set_edge_attributes(G, edge_attribs, "edge_color")
    plot = Plot(plot_width=900, plot_height=900,
                x_range=Range1d(-2.1, 2.1), y_range=Range1d(-2.1, 2.1))
    plot.title.text = "ChR2 Mouse"
    graph_renderer = from_networkx(G, nx.circular_layout, scale=1, center=(0, 0))
    pos = graph_renderer.layout_provider.graph_layout
    x, y = zip(*pos.values())
    source = ColumnDataSource({'x': pd.Series(x), 'y': pd.Series(y), 'field': list(G.nodes())})
    labels = LabelSet(x='x', y='y', text_font_size='10px', x_offset=-11, y_offset=-5, text='field', source=source)
    graph_renderer.node_renderer.glyph = Circle(size=1, fill_color=Spectral4[0])
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=0.1)
    plot.renderers.append(graph_renderer)
    plot.renderers.append(labels)
    output_file("interactive_graphs.html")
    show(plot)


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
            size= np.array(degree_list)*15,

            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
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

