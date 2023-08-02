#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:11:09 2022

@author: albit
"""


# import necessary libraries
from numpy.random import random
import networkx as nx
import pandas as pd
import numpy as np
import scipy.special as sc
import seaborn as sns
import matplotlib.patches as mpatches
import pingouin as pg
import matplotlib.pyplot as plt
#import markov_clustering as mc
from bokeh.io import output_file, show
from bokeh.models import (Circle, MultiLine, Plot, Range1d, BoxZoomTool, ResetTool)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx
from bokeh.models import ColumnDataSource, LabelSet
from IPython import get_ipython
import bokeh as bk 
import matplotlib.colors as colors
from bokeh.io import output_file, show
from bokeh.models.graphs import from_networkx
from bokeh.models import (BoxSelectTool, Circle, EdgesAndLinkedNodes, HoverTool,
                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool,)
from bokeh.models import (BoxZoomTool, Circle, HoverTool,
                          MultiLine, Plot, Range1d, ResetTool,)
from bokeh.palettes import Spectral4
#from bokeh.plotting import from_networkx
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as shc
from sklearn.decomposition import PCA
#import scola
import pylab
import networkx.algorithms.community as nx_comm
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.community import modularity
from community import community_louvain
import brainconn
from statannotations.Annotator import Annotator
import markov_clustering as mc
import random
from statsmodels.graphics.regressionplots import abline_plot
from sklearn.cluster import AgglomerativeClustering
import Data_wrangling as dw
import network_analysis.NetworkFunctions as nf
import network_analysis.plotting_utils as plu
import network_analysis.algorithms as al
import os
#%%
# simple function for loading our csv file
def loadData(data):
    data = pd.read_csv(data)  # add col list if you already made entirety of analysis usecols=(col_list)
    # data.drop(columns=data.columns[data.nunique() <= 3], inplace=True)
    data = data.apply(lambda x: x.fillna(x.mean()), axis=0)
    node_names = data.columns.to_list()
    node_number = list(item for item in range(0, len(node_names)))
    nodes = {node_number[i]: node_names[i] for i in range(len(node_number))}
    return data, nodes
#%%
#data, nodes = loadData('/Users/albitcabanmurillo/Downloads/CFOS_networks/NWX_all_avgs_n8_nonerves_c.csv')
data, nodes = loadData('/Users/albit/Desktop/Networks/CFOS_networks/NWX_all_avgs_n8_nonerves_c.csv')

bilateral = False
ROI = False
processed = False
mean = True
percentage_of_change = False
corr_plot = False
percen = np.arange(0.2, 0.8, 0.2)
percent = [round(n, 2) for n in percen]
used_P = '0.4'
sex = 'F'
#percentile = [0.2]

colors = {'KET30': 'orange', 'SAL': 'grey'}

colors = {'Cerebellum': 'grey',
 'Cortical plate': 'red',
 'Cortical subplate': 'red',
 'Hypothalamus': 'purple',
 'Medulla': 'orange',
 'Midbrain': 'green',
 'Pallidum': 'blue',
 'Pons': 'orange',
 'Striatum': 'blue',
 'Thalamus': 'violet'}

animal_dict = {'KET30': ['M1', 'M3', 'M5', 'M7', 'M9', 'M11'],
          'SAL': ['M2', 'M4', 'M6', 'M8', 'M10', 'M12']}

filepath = 'C:/Users/cds4/Desktop/Network'
Allen_areas_dict_filepath = "C:/Users/cds4/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle"
results_folder = 'C:/Users/cds4/Desktop/Network/LabMeeting'
data_wrangled, filepath, results_folder, data_mean = dw.dataManager(filepath, results_folder, bilateral, ROI, processed, mean)
#extra_discard = ['MB', 'MY', 'P']
#for e in data_wrangled.keys():
 #   data_wrangled[e] = data_wrangled[e].drop(extra_discard, axis = 1)
if ROI:
    data_wrangled = {'KET30':  'C:/Users/cds4/Desktop/Network/KET30_Large_Network_ROI.csv',
                 'SAL':'C:/Users/cds4/Desktop/Network/SAL_Large_Network_ROI.csv' }
else:
    data_wrangled = {'KET30':  'C:/Users/cds4/Desktop/Network/KET30_Large_Network.csv',
                 'SAL':'C:/Users/cds4/Desktop/Network/SAL_Large_Network.csv' }
ROIs, Allen_Groups = dw.ag_pickle(Allen_areas_dict_filepath)
ROIs_dict = dw.get_ROIs_dict(filepath)
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
nodes = {exp: nf.loadData(data_wrangled[exp], csv = True)[1] for exp in data_wrangled.keys()}


#%%
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
#%%
corr, p = corrMatrix(data)

#%%
correlation_matrix ={exp: nf.corrMatrix(raw_data[exp], z_trans=True) for exp in raw_data.keys()}



# using this function we will threshold based off of p-values previously calculated
def significanceCheck(p, corr, alpha, threshold=0.0, names=None, title = '', results_folder= '',plot=True, include_Negs=True):
    p = np.where((p >= alpha), 0, p)  # if not significant --> zero
    p = np.where((p != 0), 1, p)  # if significant --> one
    if not include_Negs:
        p = np.where((corr < 0), 0, p)
    threshold_matrix = np.multiply(corr, p)  # remove any insignificant correlations
    # remove correlations below threshold
    threshold_matrix = np.where((abs(threshold_matrix) < threshold), 0, threshold_matrix)
    # create a heatmap of correlations if wanted
    if plot:
        pandas_matrix = pd.DataFrame(threshold_matrix, index=list(names.values()), columns=list(names.values()))
        sns.set(font_scale=2.5)
        sns.clustermap(pandas_matrix, cmap="rocket", figsize=(35, 35), cbar_pos=(.02,.85,.05,.18))
        plt.savefig(os.path.join(results_folder, 'Euclidean_corr_matrix_' + title + '.png'))
    return threshold_matrix, pandas_matrix 

#%%

sig, sigData = significanceCheck(p, corr, 0.01, names = nodes, threshold = 0.5, include_Negs=(False))

threshold_matrix = {exp: significanceCheck(matrix[2], matrix[0],
                                              alpha = 1, names = nodes[exp], title = exp, results_folder=results_folder, include_Negs = False)\
                    for exp, matrix in correlation_matrix.items()}
# Function that will threshold R-vals on the top percentile
def percentile(array, p, names = None, title = '', results_folder= '', plot = True, Anatomy = None):
    num_obs = int(np.size(array, 0) ** 2 * p)
    crit_value = -np.sort(-array.flatten())[num_obs - 1]
    percent_arr = np.where(array < crit_value, 0, array)
    pandas_matrix = pd.DataFrame(percent_arr, index=list(names.values()), columns=list(names.values()))
    if plot:
        sns.set(font_scale=2.5)
        sns.clustermap(pandas_matrix, cmap="rocket", figsize=(35, 35), cbar_pos=(.02,.8,.05,.18))
        plt.savefig(os.path.join(results_folder, 'Euclidean_corr_matrix_' + title + '.png'))
        if Anatomy:
            sorted_dict = dict(sorted(Anatomy.items(), key=lambda item: item[1]))
            list_for_sort = list(sorted_dict.keys())

            # Reorganize the pandas_matrix to reflect the order of the Allen ROIs
            allen_pandas = pandas_matrix[list_for_sort].loc[list_for_sort]

            # Create color assignments for Allen ROIs
            num_allens = list(sorted_dict.values())
            allens_unique = np.unique(num_allens)
            color_list = [color for color in sns.color_palette('Paired', n_colors=len(allens_unique))]
            color_ref = dict(zip(map(str, allens_unique), color_list))
            allen_colors = pd.Series(num_allens, index=allen_pandas.columns).map(color_ref)

            # Create a legend for the Allen ROIs
            cerebellum = mpatches.Patch(color=color_list[0], label='Cerebellum')
            cort_plate = mpatches.Patch(color=color_list[1], label='Cortical Plate')
            cort_subplate = mpatches.Patch(color=color_list[2], label='Cortical Subplate')
            hypothalamus = mpatches.Patch(color=color_list[3], label='Hypothalamus')
            medulla = mpatches.Patch(color=color_list[4], label='Medulla')
            midbrain = mpatches.Patch(color=color_list[5], label='Midbrain')
            pallidum = mpatches.Patch(color=color_list[6], label='Pallidum')
            pons = mpatches.Patch(color=color_list[7], label='Pons')
            striatum = mpatches.Patch(color=color_list[8], label='Striatum')
            thalamus = mpatches.Patch(color=color_list[9], label='Thalamus')

            # Plot the newly generated Allen ROI correlation maitrx
            plt.figure()
            sns.clustermap(allen_pandas, cmap='rocket', row_colors=allen_colors, col_colors=allen_colors,
                           row_cluster=False, col_cluster=False, xticklabels=False, yticklabels=False,
                           figsize=(35, 35), cbar_pos=(.02,.8,.05,.18), cbar_kws={'label': 'arctan(R)'})
            plt.legend(
                handles=[cerebellum, cort_plate, cort_subplate, hypothalamus, medulla, midbrain, pallidum, pons,
                         striatum, thalamus],
                bbox_to_anchor=(19.0, 1.13), prop = {'size': 35})
            plt.savefig(os.path.join(results_folder, 'Anatomical_corr_matrix_' + title + '.png'))
            plt.show()
    return percent_arr, pandas_matrix

per_adj = {exp: {str(p): percentile(matrix[0], p, names = nodes[exp], title = '_'.join([exp, str(p)]), results_folder = results_folder, plot = True,Anatomy= ROIs) for p in percent}\
                    for exp, matrix in threshold_matrix.items()}

 #%%
# we will create our undirected network graphs based on our matrices
def networx(corr_data, nodeLabel):
    graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)
    graph = nx.relabel_nodes(graph, nodeLabel)
    pos = nx.circular_layout(graph)
    nx.set_node_attributes(graph, pos, name='pos')
    return graph, pos
#%%

G,y = networx(sig, nodes)
graph = {exp: {p: nf.networx(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}

#%%

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
    labels = LabelSet(x='x', y='y', text_font_size='10px', x_offset=0, y_offset=0, text='field', source=source)
    graph_renderer.node_renderer.glyph = Circle(size=50, fill_color=Spectral4[0])
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=1)
    plot.renderers.append(graph_renderer)
    plot.renderers.append(labels)
    output_file("interactive_graphs.html")
    show(plot)
    
#%%

GraphingNetwork(G, ".05", nx.circular_layout)

[[GraphingNetwork(G, '_'.join([exp,p]), nx.circular_layout)for p, G in pe.items()] for exp, pe in graph.items()]
#%%
sigData_2 = np.asmatrix(sigData)
sigData_2 = {exp: {p: np.asmatrix(matrix[1]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}

ROIIDdata = sigData.index.values
ROIIDdata = {exp: {p: matrix[1].index.values for p, matrix in pe.items()} for exp, pe in per_adj.items()}

sigData_3 = np.asmatrix(sig)
sigData_3 = {exp: {p: np.asmatrix(matrix[0]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}


Gd = nx.to_networkx_graph(sigData_3)    #Added this on 6/7 because I think it will remove the self edges which are messing with the clustering coeficient
Gd = {exp: {p: nx.to_networkx_graph(matrix) for p, matrix in pe.items()} for exp, pe in sigData_3.items()}

#Gd.remove_edges_from(nx.selfloop_edges(Gd))
#Gd = {exp: {p: matrix.remove_edges_from(nx.selfloop_edges(matrix)) for p, matrix in pe.items()} for exp, pe in Gd.items()}

for exp,pe in Gd.items():
    for p, matrix in pe.items():
        Gd[exp][p].remove_edges_from(nx.selfloop_edges(Gd[exp][p]))

#relabels the nodes to match the cell names
Gd = nx.relabel_nodes(Gd,lambda x: ROIIDdata[x])
Gd_new = {exp: {p: nx.relabel_nodes(matrix,lambda x: ROIIDdata[exp][p][x]) for p, matrix in pe.items()} for exp, pe in Gd.items()}

#%%

#for gephi use create G through this function instead

def create_corr_network_plz(G, title = '', results_folder = ''):
    #crates a list for edges and for the weights
    
    edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())


   # weights = tuple([(5*(x))**2 for x in weights]) 
    
    #positions, how to make the graph shpae
    positions=nx.circular_layout(G)
    
    #Figure size
    plt.figure(figsize=(15,15))

    #draws nodes
    nx.draw_networkx_nodes(G,positions,node_color='b',
                           node_size=500,alpha=0.8)
    
    #Styling for labels
    nx.draw_networkx_labels(G, positions, font_size=8, 
                            font_family='sans-serif')
        
    #draws the edges
    nx.draw_networkx_edges(G, positions, edgelist=edges,style='solid', width = 10*(weights), )
    
    # displays the graph without axis
    plt.axis('off')
    #saves image
    plt.savefig(os.path.join(results_folder, 'network' + title + '.png'))
    #plt.savefig("part1.png", format="PNG")
    plt.show()


#%%
create_corr_network_plz(Gd)

[[create_corr_network_plz(matrix, title = '_'.join([exp, str(p)]), results_folder = results_folder) for p, matrix in pe.items()] for exp, pe in Gd.items()]

def gephiMyNetwork(sig, exp, p):
    Gg = nx.to_networkx_graph(sig)  #I dont know why but gephi function only works when we make the network through here, but everything is the same.
 #   path = ‘/Users/albitcabanmurillo/Downloads/CFOS_networks/Networks_clean/Gephi_a01r07_fullnetwork.gexf’
    path = os.path.join(results_folder, 'Participation_coeficient_vs_Module_Degree' + exp  + p +'.gexf')
    nx.write_gexf(Gg, path)
 #   resultstest.to_csv(r’/Users/albitcabanmurillo/Downloads/CFOS_networks/Networks_clean//Hubscores_test.csv’, index = True, header=True)
#%%
for exp, s in per_adj.items():
    for p, g in s.items():
        gephiMyNetwork(g[0], exp, p)

#%%
#calc shortest path cause for some reason its having trouble with it inside the function
def shortest(G, sigData):

    short = nx.floyd_warshall_numpy(G, weight='weight')
    shortavg = np.mean(short, axis = 0)
    keys = sigData.index.values.tolist()
    vals = shortavg.tolist()
    zip_iterator = zip(keys, vals)
    short_dictionary = dict(zip_iterator)
    
    return short_dictionary


#%%

short_avg = shortest(G)

short_avg = {exp: {p: shortest(matrix, per_adj[exp][p][1]) for p, matrix in pe.items()} for exp, pe in graph.items()}

#%%

def FindMyHubs(G, sigData):
    
    G_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G.edges(data='weight')}  # creates a dict of calculted distance between all nodes
    nx.set_edge_attributes(G, values= G_distance_dict, name='distance') #sets the distance as an atribute to all nodes
    
    cluster_coefficient = nx.clustering(G, weight='weight')     #calculated different measures of importance
    degree_cent = nx.degree_centrality(G)
    eigen_cent = nx.eigenvector_centrality(G, max_iter=100000, weight='weight') 
    betweenness = nx.betweenness_centrality(G, weight='distance', normalized=True)
    closeness = nx.closeness_centrality(G, distance='distance')
    degree = list(G.degree(G, weight='weight'))
    communicability = nx.communicability_betweenness_centrality(G)
    short_avg = shortest(G, sigData)

    dict = {'eigen_cent': eigen_cent, 'betweenness': betweenness, 'closeness': closeness, 
            'clustering':cluster_coefficient, 'degree_cent':degree_cent, 'communicability': communicability, 'short_avg':short_avg} 

    Results = pd.DataFrame(dict)    #create a frame with all the measures of importance for every region
    
    #Van den Huevel(2010) - https://www.jneurosci.org/content/30/47/15915
    #used the top or bottom quartiles to determine the hubness of all nodes so here we calculate that.
    
    eig_quart_upper = Results.eigen_cent.quantile(0.75)  #eig centrality upper and lower quartiles
    #eig_quart_lower = Results.eigen_cent.quantile(.25)
    degree_cent_quart_upper = Results.degree_cent.quantile(0.75) #weighted degree centrality upper quartiles
    betweeness_quart_upper = Results.betweenness.quantile(0.75) #betweeness centrality upper quartile 
    closeness_quart_upper = Results.closeness.quantile(.75)    #closeness centrality upper quartile
    clustering_quart_lower = Results.clustering.quantile(.25)  #clustering lower quartile
    com_quart_upper = Results.communicability.quantile(.75)    #communicability upper quartile
    short_quart_lower = Results.short_avg.quantile(.25)
    quartiles_dict = {'eing_cent': eig_quart_upper, 'betweeness': betweeness_quart_upper, 'clustering': clustering_quart_lower, 'shortness':  short_quart_lower}
    #for each significant measure an ROI has add one in the score column, a score >= 2 is considered a hub node.
    Results['score'] = 0
    Results['score'] = np.where((Results['eigen_cent'] > eig_quart_upper), Results['score']+1, Results['score'])
    #Results['score'] = np.where((Results['eigen_cent'] < eig_quart_lower), Results['score']+1, Results['score'])
    Results['score'] = np.where((Results['betweenness'] >= betweeness_quart_upper), Results['score']+1, Results['score'])
    Results['score'] = np.where((Results['clustering'] <= clustering_quart_lower), Results['score']+1, Results['score'])
    #Results['score'] = np.where((Results['communicability'] >= com_quart_upper), Results['score']+1, Results['score'])
    Results['score'] = np.where((Results['short_avg'] <= short_quart_lower), Results['score']+1, Results['score'])
    
    NonHubs =  Results[(Results['score'] < 2) ].index   #create an index of rois with a score of less than 2 in hubness
    
    Hubs = Results.drop(NonHubs, errors = 'ignore')    #create a new frame with only the important nodes
    
    return Results, Hubs, quartiles_dict
#%%

resultstest, resultsbsidetest = FindMyHubs(G) #hubs are bside

results = {exp: {p: FindMyHubs(matrix , per_adj[exp][p][1]) for p, matrix in pe.items()} for exp, pe in graph.items()}
dw.pickle_save(results, results_folder + '/results.pickle')
KetHubs = np.histogram(np.array(results['KET30']['0.4'][0].score), np.arange(0,6,1))
SalHubs = np.histogram(np.array(results['SAL']['0.4'][0].score), np.arange(0,6,1))

fig = plt.figure(figsize=(14,10))
fig.patch.set_facecolor('white')
plt.plot(np.arange(5), KetHubs[0], color = 'orange', marker = 'o', linewidth = 5, markeredgewidth = 8)
plt.plot(np.arange(5), SalHubs[0], color = 'grey', marker = 'o', linewidth = 5, markeredgewidth = 8)
plt.vlines(x = 1.5, ymin = 0, ymax = 80, color = 'r')
plt.xticks([0,1,2,3,4], ['0','1','2','3','4'], fontsize =30)
plt.yticks([0,20,40,60], ['0','20','40','60'], fontsize =30)
plt.xlabel('Score', fontsize =30)
plt.ylabel('Hubs', fontsize =30)

fig.legend(handles = plu.create_custom_legend({'KET30': 'orange', 'SAL':'grey'}), prop={"size":25})
#plt.legend(loc='upper right')
plt.show()
fig.savefig(os.path.join(results_folder, 'Number of hubs' + exp + '_' + p + '.png'), dpi=300)

for exp, pe in results.items():
    for p, matrix in pe.items():
        fig, ax = plt.subplots(nrows=1, ncols=4, gridspec_kw={
            'hspace': 0.4, 'wspace': 0.43}, figsize=(12,8))
        sns.violinplot(y = 'eigen_cent', data = matrix[0], ax = ax[0], inner = 'points')
        ax[0].hlines(y = matrix[2]['eing_cent'], xmin = ax[0].get_xlim()[0], xmax =ax[0].get_xlim()[1], color = 'red')
        ax[0].set_title('Eigen centrality', fontsize =20)
        ax[0].set_ylim(-0.05, 0.25)
        ax[0].set_yticks([0, 0.10, 0.20], ['0', '0.10', '0.20'], fontsize = 15)
        sns.violinplot(y = 'betweenness', data = matrix[0], ax = ax[1], inner = 'points', color = 'orange')
        ax[1].hlines(y = matrix[2]['betweeness'], xmin = ax[1].get_xlim()[0], xmax =ax[1].get_xlim()[1], color = 'red')
        ax[1].set_title('Betweenness', fontsize =20)
        ax[1].set_ylim(-0.005, 0.035)
        ax[1].set_yticks([0, 0.015, 0.03], ['0', '0.015', '0.035'], fontsize = 15)
        sns.violinplot(y = 'clustering', data = matrix[0], ax = ax[2], inner = 'points', color = 'purple')
        ax[2].hlines(y = matrix[2]['clustering'], xmin = ax[2].get_xlim()[0], xmax =ax[2].get_xlim()[1], color = 'red')
        ax[2].set_title('Clustering', fontsize =20)
        ax[2].set_ylim(-0.05, 0.4)
        ax[2].set_yticks([0, 0.15, 0.30], ['0', '0.15', '0.30'], fontsize = 15)
        sns.violinplot(y = 'short_avg', data = matrix[0], ax = ax[3], inner = 'points', color = 'green')
        ax[3].hlines(y = matrix[2]['shortness'], xmin = ax[3].get_xlim()[0], xmax =ax[3].get_xlim()[1], color = 'red')
        ax[3].set_title('Shortness', fontsize =20)
        ax[3].set_ylim(0.45, 2.2)
        ax[3].set_yticks([0.5, 1, 2], ['0.5', '1', '2'], fontsize = 15)
        ax[0].set_ylabel('')
        ax[1].set_ylabel('')
        ax[2].set_ylabel('')
        ax[3].set_ylabel('')
        fig.suptitle('Stadistical analysis Hubs ' + exp + ' ' + p )
        plt.show()
        fig.savefig(os.path.join(results_folder, 'Stadistical_analysis_Hubs_'+ exp + '_' + p  + '.png'), dpi=300)

#%%
#Use your clustering of preference here I just added a simple louvain one to make the point just 
#make it into an array for everything else

from community import community_louvain

louv = community_louvain.best_partition(G)

louv = {exp: {p: community_louvain.best_partition(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}
#%%

Modules_full_network = np.array(list(louv.items()))
Modules= {exp: {p: np.array(list(matrix.items()))[:,1] for p, matrix in pe.items()} for exp, pe in louv.items()}
Modules = Modules_full_network[:,1]

#%% function to calculate module statistics

#%%
indexx = resultsbsidetest.index.copy()
index = {exp: {p: matrix[1].index.copy() for p, matrix in pe.items()} for exp, pe in results.items()}
indexx_full_network = resultstest.index.copy()
indexx_full_network = {exp: {p: matrix[0].index.copy() for p, matrix in pe.items()} for exp, pe in results.items()}

#modules = pd.DataFrame(Modules, columns=['Module'],index = indexx )
modules_full_network = pd.DataFrame(Modules, columns=['Module'],index = indexx_full_network )

modules_full_network = {exp: {p: pd.DataFrame(matrix, columns=['Module'],index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in Modules.items()}


#within_module_degree = brainconn.centrality.module_degree_zscore(sigc, modules)
within_module_degree_FN = brainconn.centrality.module_degree_zscore(sig, modules_full_network)

within_module_degree_FN = {exp: {p:  brainconn.centrality.module_degree_zscore(matrix[0], modules_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}


#Ppos, Pneg = brainconn.centrality.participation_coef_sign(sigc, modules)
Ppos_FN, Pneg_FN = brainconn.centrality.participation_coef_sign(sig, modules_full_network)

FN = {exp: {p:  brainconn.centrality.participation_coef_sign(matrix[0], modules_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}


#Z_data = pd.DataFrame(within_module_degree, columns=['Within module degree'],index = indexx )
#P_data = pd.DataFrame(Ppos, columns=['Participation coefficient'], index = indexx)
  
Z_data_FN = pd.DataFrame(within_module_degree_FN, columns=['Within module degree'],index = indexx_full_network )
P_data_FN = pd.DataFrame(Ppos_FN, columns=['Participation coefficient'], index = indexx_full_network)

Z_data_FN = {exp: {p:  pd.DataFrame(matrix, columns=['Within module degree'],index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in within_module_degree_FN.items()}
P_data_FN = {exp: {p:  pd.DataFrame(matrix[0], columns=['Participation coefficient'], index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in FN.items()}

pc_df = pd.read_csv("C:/Users/cds4/Desktop/Network/\RESULTS_COMBINED/pc_dfKET30.csv")
pc_df = pc_df.set_index('Unnamed: 0')


#Hub_all_data = pd.concat([resultsbsidetest, Z_data], axis=1 )
#Hub_all_data = pd.concat([Hub_all_data, P_data], axis = 1)
#Hub_all_data = pd.concat([Hub_all_data, modules], axis = 1)

Network_all_data = pd.concat([resultstest, Z_data_FN], axis=1 )
Network_all_data = pd.concat([Network_all_data, P_data_FN], axis = 1)
Network_all_data = pd.concat([Network_all_data, modules_full_network], axis = 1)

Network_all_data = {exp: {p: pd.concat([matrix[0], Z_data_FN[exp][p], P_data_FN[exp][p], modules_full_network[exp][p], pc_df], axis=1) for p, matrix in pe.items()} for exp, pe in results.items()}


Hubs_per_area = {}
Network_area_combined = {}
for exp, pe in Network_all_data.items():
    for p,m in pe.items():

        Network_all_data[exp][p] = m[['Complete name', 'Region', 'score', 'Module', 'Percentage of change', 'p-value', 'Within module degree',
               'Participation coefficient','eigen_cent', 'betweenness', 'closeness', 'clustering', 'degree_cent',
               'communicability', 'short_avg']]

        df = pd.DataFrame.from_dict({r: [sum(m.query('Region ==@r').score>1)] for r in sorted(set(m.Region))})
        df['exp'] = exp
        ndata = Network_all_data[exp][p].copy()
        ndata['exp'] = exp
        if p not in Hubs_per_area.keys():
            Hubs_per_area[p] = []
            Network_area_combined[p] = []
        Hubs_per_area[p].append(df)
        Network_area_combined[p].append(ndata)

for p in Hubs_per_area.keys():
    Hubs_per_area[p] = pd.concat(Hubs_per_area[p])
    Network_area_combined[p] = pd.concat(Network_area_combined[p])

data_stats = Network_area_combined[used_P]
stats = []
no_stats = ['Complete name', 'Region', 'score', 'Module', 'Percentage of change',
       'p-value', 'exp']

for col in list(data_stats.columns):
    if col not in no_stats:
        n = np.array(data_stats[col])
        normality = pg.normality(n)
        param = normality['normal'][0]
        ket = data_stats.query('exp == "KET30"')
        sal = data_stats.query('exp == "SAL"')
        x = np.array(ket[col])
        y = np.array(sal[col])
        if param:
            df = pg.ttest(x,y)
            df['Measure'] = col
            stats.append(df)
        else:
            df = pg.mwu(x, y)
            df['Measure'] = col
            stats.append(df)
stats = pd.concat(stats)


data_toplot = Hubs_per_area[used_P]
fig = plt.figure(figsize=(14,10))
fig.patch.set_facecolor('white')
rows = data_toplot.loc[0, :].values.tolist()
for row in rows:
    exp = row[-1]
    plt.plot(np.arange(len(row) - 1),row[:-1], color=colors_exp[exp], marker='o', linewidth=5, markeredgewidth=8)
    plt.xticks(np.arange(len(row) - 1), list(data_toplot.columns[:-1]), fontsize=30, rotation = 90)
  #  plt.yticks([0, 20, 40, 60], ['0', '20', '40', '60'], fontsize=30)
    plt.xlabel('Regions', fontsize=30)
    plt.ylabel('Hubs', fontsize=30)


fig.legend(handles=plu.create_custom_legend({'KET30': 'orange', 'SAL': 'grey'}), prop={"size": 25}, loc = 4)
# plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
fig.savefig(os.path.join(results_folder, 'Number of hubs per area' + exp + '_' + used_P + '.png'), dpi=300)

[[df.to_csv(os.path.join(results_folder, 'Score_All' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Network_all_data.items()]

[[df.rename(columns = {'Module': 'Modularity'}).to_csv(os.path.join(results_folder, 'gephi' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Network_all_data.items()]


shared_hubs = []
hubs_exp = []
for exp,per in Network_all_data.items():
    hubs_exp.append(list(per[used_P][per[used_P].score > 1].index))
shared_hubs = set(hubs_exp[0]).intersection(hubs_exp[1])

#%%
#plot z vs p in pretty scatter plot but not bu clusters, only shows all data in one color.
for exp, pe in Network_all_data.items():
    for p, matrix in pe.items():
        fig = plt.figure(figsize=(20, 20))
        Modules = matrix['Module'].to_numpy()
        my_col = {1: "lightcoral", 2: "lightskyblue", 3: "lightgreen"}
        fig = sns.JointGrid(data=matrix, y = "Within module degree", x= "Participation coefficient")
        fig.plot_joint(sns.scatterplot, s=1, color='cornflowerblue')
        fig.plot_marginals(sns.kdeplot, color='cornflowerblue', fill=True)
        fig.ax_marg_x.set_xlim(ax.get_xlim()[0], 0.8) #calcular percentiles y cambiar
        fig.ax_marg_y.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1])
        fig.set_axis_labels('P', 'Z', fontsize=22)
        ax.vlines(.3, ymin = ax.get_ylim()[0], ymax= ax.get_ylim()[1],color="black", linestyle="--")

        for idx, (x,y) in matrix.loc[:, ["Participation coefficient" , "Within module degree"]].iterrows():
            fig.ax_joint.annotate(idx, (x,y), fontsize=9, weight = 'bold')

        for q in np.quantile(matrix['Participation coefficient'], [0.25, 0.75]):
             for ax in (fig.ax_joint, fig.ax_marg_x):
                  ax.axvline(q, color="black", linestyle="--")
        for q in np.quantile(matrix['Within module degree'], [0.25, 0.75]):
             for ax in (fig.ax_joint, fig.ax_marg_y):
                  ax.axhline(q, color="black", linestyle="--")
     #   ax.set_title('Participation coeficient vs Module Degree' + exp + ' ' + p)
       # fig.suptitle('Participation coeficient vs Module Degree' + exp + ' ' + p)
        fig.savefig(os.path.join(results_folder, 'Participation_coeficient_vs_Module_Degree' + exp + '_' + p +'.png'), dpi=300)


#%%
#plot z vs p in pretty scatter plot but not bu clusters, only shows all data in one color.
fig = plt.figure(figsize=(20, 20))
Modules = Network_all_data['Module'].to_numpy()
my_col = {1: "lightcoral", 2: "lightskyblue", 3: "lightgreen"}
g = sns.JointGrid(data=Network_all_data, y = "Within module degree", x= "Participation coefficient")
g.plot_joint(sns.scatterplot, s=50, color='cornflowerblue')
g.plot_marginals(sns.kdeplot, color='cornflowerblue', fill=True)
#g.ax_marg_x.set_xlim(.3, .85)
#g.ax_marg_y.set_ylim(.70, 2)
g.set_axis_labels('P', 'Z', fontsize=20)

#Plot the names of the rois 

for idx, (x,y) in Network_all_data.loc[:, ["Participation coefficient" , "Within module degree"]].iterrows():
    g.ax_joint.annotate(idx, (x,y), fontsize=8, weight = 'bold') 


for ax in (g.ax_joint, g.ax_marg_x):
   ax.axvline(.3, color="black", linestyle="--")

for ax in (g.ax_joint, g.ax_marg_x):
   ax.axvline(.75, color="black", linestyle="--", )

#for ax in (g.ax_joint, g.ax_marg_x):
 #  ax.axvline(.75, color="black", linestyle="--")

for q in np.quantile(Network_all_data['Within module degree'], [0.80, 0.20]):
     for ax in (g.ax_joint, g.ax_marg_y):
          ax.axhline(q, color="black", linestyle="--")
          

similarity = [al.Similarity(per_adj['KET30'][str(p)][0], per_adj['SAL'][str(p)][0]) for p in percent]
