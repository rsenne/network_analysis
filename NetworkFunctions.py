"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by cesar coelho and gisella vetere.
we thank them for their kind support throughout this process"""
# author:ryan senne/ramirez group

# import necessary libraries
import networkx as nx
import pandas as pd
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

file = '/home/ryan/The_Ramirez_Group/Cell_Counts.csv'


# simple function for loading our csv file
def loadData(file):
    data = pd.read_csv(file)
    return data


# correlate our c-Fos counts between brain regions, df for data
# type for correlation coefficient i.e. "pearson"
def corrMatrix(data):
    r = np.corrcoef(data, rowvar=False)
    rf = r[np.triu_indices(r.shape[0], 1)]
    df = data.shape[1] - 2
    ts = rf * rf * (df / (1 - rf * rf))
    pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))
    p = np.zeros(shape=r.shape)
    p[np.triu_indices(p.shape[0], 1)] = pf
    p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]
    p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
    return r, p


# using this function we will threshold based off of p-values previously calculated
def significanceCheck(p, corr, alpha):
    p = np.where((corr >= alpha), 0, p)
    p = np.where((p != 0), 1, p)
    thresholded_matrix = np.multiply(corr, p)
    return thresholded_matrix


# we will create our undirected network graphs based on our matrices
def networx(corr_data):
    graph = nx.from_numpy_array(corr_data)
    plt.subplot()
    nx.draw(graph, with_labels=True, node_color="Blue", node_size=200, edge_color="black", linewidths=2, font_size=15)
    plt.show()
    return graph


# calculate the necessary graph attributes such as centrality, betweenness, global efficiency etc.
def graphAttributes(graph):
