"""this is a set of functions necessary for the creation of undirected c-Fos networks.
this project was inspired and adapted from work done by Drs. Cesar Coelho, Anne Wheeler, and Gisella Vetere.
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
import matplotlib.patches as mpatches
from bct.algorithms import centrality


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
            pandas_matrix = pd.DataFrame(threshold_matrix,index=list(names.values()),columns=list(names.values()))
            if Anatomy:
                #Create a sorted dictionary from the unpickled ROIs dictionary
                sorted_dict = dict(sorted(Anatomy.items(),key=lambda item: item[1]))
                list_for_sort = list(sorted_dict.keys())

                #Reorganize the pandas_matrix to reflect the order of the Allen ROIs
                allen_pandas = pandas_matrix[list_for_sort].loc[list_for_sort]

                #Create color assignments for Allen ROIs
                num_allens = list(sorted_dict.values())
                allens_unique = np.unique(num_allens)
                color_list = [color for color in sns.color_palette('Set3',n_colors=len(allens_unique))]
                color_ref = dict(zip(map(str,allens_unique),color_list))
                allen_colors = pd.Series(num_allens,index=allen_pandas.columns).map(color_ref)

                #Create a legend for the Allen ROIs
                cerebellum = mpatches.Patch(color=color_list[0],label='Cerebellum')
                cort_plate = mpatches.Patch(color=color_list[1],label='Cortical Plate')
                cort_subplate = mpatches.Patch(color=color_list[2],label='Cortical Subplate')
                hypothalamus = mpatches.Patch(color=color_list[3],label='Hypothalamus')
                medulla = mpatches.Patch(color=color_list[4],label='Medulla')
                midbrain = mpatches.Patch(color=color_list[5],label='Midbrain')
                pallidum = mpatches.Patch(color=color_list[6],label='Pallidum')
                pons = mpatches.Patch(color=color_list[7],label='Pons')
                striatum = mpatches.Patch(color=color_list[8],label='Striatum')
                thalamus = mpatches.Patch(color=color_list[9],label='Thalamus')

                #Plot the newly generated Allen ROI correlation maitrx
                plt.figure()
                sns.clustermap(allen_pandas,cmap='viridis',row_colors=allen_colors,col_colors=allen_colors,
                               row_cluster=False,col_cluster=False,xticklabels=False,yticklabels=False,
                               figsize=(10,10),cbar_pos=(0.1,0.15,.02,.4),cbar_kws={'label':'Pearson Correlation (R)'})
                plt.legend(handles=[cerebellum,cort_plate,cort_subplate,hypothalamus,medulla,midbrain,pallidum,pons,striatum,thalamus],
                           bbox_to_anchor=(5.0,1.6))
        else:
            pandas_matrix = pd.DataFrame(threshold_matrix)
        sns.clustermap(pandas_matrix,cmap='viridis',method='ward',metric='euclidean',figsize=(10,10),cbar_pos=(.9,.9,.02,.10))
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


def cluster_attributes(graph,communities):
    adj_matrix = nx.to_numpy_matrix(graph) #Will create an adjacency matrix from the graph
    WMDz = centrality.module_degree_zscore(adj_matrix,communities,flag=0) #calculate the WMDz
    PC = centrality.participation_coef(adj_matrix,communities,'undirected')
    return WMDz,PC


def get_ordered_degree_list(G):
    degree_ordered = {k: v for k, v in sorted(dict(G.degree()).items(), key=lambda item: item[1])}
    return degree_ordered
