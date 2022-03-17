"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

# import necessary libraries

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.special as sc
import seaborn as sns
from statsmodels.sandbox.stats.multicomp import multipletests


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


def grab_attributes(graph):
    deg = nx.degree_centrality(graph)
    between = nx.betweenness_centrality(graph)
    eig = nx.eigenvector_centrality(graph)
    deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[1])}
    between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[1])}
    eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[1])}
    return deg_sort, between_sort, eig_sort


def get_ordered_degree_list(G):
    degree_ordered = {k: v for k, v in sorted(dict(G.degree()).items(), key=lambda item: item[1])}
    return degree_ordered
