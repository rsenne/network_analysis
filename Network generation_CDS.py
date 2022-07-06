import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import Data_wrangling as dw
import network_analysis.NetworkFunctions as nf
import network_analysis.plotting_utils as plu
import network_analysis.algorithms as al
import Quantification as q
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import pandas as pd
bilateral = False
ROI = False
processed = False
mean = True
percentage_of_change = False
corr_plot = False
percentile = np.arange(0.05, 1, 0.05)
percentile = [round(n, 2) for n in percentile]
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


filepath = 'C:/Users/cds4/Desktop/Network'
Allen_areas_dict_filepath = "C:/Users/cds4/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle"
results_folder = 'C:/Users/cds4/Desktop/Network/Results'
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


if percentage_of_change:
    #Areas ordered anatomically
    norm_group = 'SAL'
    order_areas = data_mean['SAL'].columns
    order_complete_name = [ROIs_dict[ab] for ab in order_areas]
    order_br = [ROIs[a] for a in order_areas]
    colors = plu.color_allens(order_br)
    color_order = [colors[b] for b in order_br]
    #percentage of change with respect of the baseline
    percentage_change_dict = q.percentage_change(data_mean, norm_group = 'SAL')
    stats_dict = {e: q.tstud(data, data_wrangled[norm_group], threshold = 0.05) for e, data in data_wrangled.items() if e != norm_group}
    pc_df = {}
    sig_pc_df = {}
    sig_pc_dict = {}
    order_sig_br = {}
    for e in percentage_change_dict.keys():
        pc_df[e] = pd.concat([percentage_change_dict[e][2], stats_dict[e][0]], axis=1)
        cn = []
        reg = []
        sig_df = []
        for r in pc_df[e].iterrows():
            cn.append(ROIs_dict[r[0]])
            reg.append(ROIs[r[0]])
        pc_df[e]['Complete name'] = cn
        pc_df[e]['Region'] = reg
        pc_df[e] = pc_df[e][['Complete name', 'Region','Percentage of change', 'p-value']]
        pc_df[e] = pc_df[e].sort_values('Region')
        sig_df = pc_df[e].loc[stats_dict[e][1]]
        sig_values = np.array(sig_df['Percentage of change'])
        reg_str = np.array(sig_df['Region'])
        colors_sig = [colors[r] for r in list(sig_df['Region'])]
        order_sig_br[e] = [np.arange(len(list(sig_df['Region']))),list(sig_df['Complete name'])]

        sig_pc_df[e] = [sig_values, reg_str, colors_sig, sig_df]

    [sig_pc_df[e][3].to_csv(os.path.join(results_folder, 'sig_pc_df' + e + \
                                 '.csv')) for e in sig_pc_df.keys()]
    [pc_df[e].to_csv(os.path.join(results_folder, 'pc_df' + e + \
                                 '.csv')) for e in pc_df.keys()]

    x_pos = np.arange(len(order_br))

    fig,ax = plt.subplots(figsize=(100,50))
    for e in percentage_change_dict.keys():
        ax.bar(x_pos, percentage_change_dict[e][0], width = 0.6, color= color_order, tick_label = percentage_change_dict[e][1])
    ax.set_xticklabels(percentage_change_dict[e][1], rotation = 90, fontsize= 30)
    #ax.set_xticklabels(order_complete_name, rotation = 90, fontsize= 30)
    ylabels = ax.get_yticks().tolist()
    label_format = '{:,.0f}'
    ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
    ax.set_yticklabels([label_format.format(x) for x in ylabels],  fontsize= 60)
    fig.legend(handles = plu.create_custom_legend(colors),  fontsize = 'xxx-large', prop={"size":50})
    title = 'Percentage of change' + e
    #ax.set_title(title, fontsize = 70)
    ax.set_ylabel('Percentage of change', fontsize = 70)
    fig.savefig(os.path.join(results_folder, title+'.png'))


    fig,ax = plt.subplots(figsize=(75,50))
    for e in sig_pc_df.keys():
        ax.bar(order_sig_br[e][0], sig_pc_df[e][0], width = 0.6, color= sig_pc_df[e][2], tick_label = sig_pc_df[e][1])
    #ax.set_xticklabels(percentage_change_dict[e][1], rotation = 90, fontsize= 30)
    ax.set_xticklabels(order_sig_br[e][1], rotation = 90, fontsize= 30)
    ylabels = ax.get_yticks().tolist()
    label_format = '{:,.0f}'
    ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
    ax.set_yticklabels([label_format.format(x) for x in ylabels],  fontsize= 60)
    fig.legend(handles = plu.create_custom_legend(colors),  fontsize = 'xx-large', prop={"size":40})
    title = 'Percentage of significant change ' + e
    #ax.set_title(title, fontsize = 70)
    ax.set_ylabel('Percentage of change', fontsize = 70)
    fig.savefig(os.path.join(results_folder, title+'.png'))



#Network generation
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
nodes = {exp: nf.loadData(data_wrangled[exp], csv = True)[1] for exp in data_wrangled.keys()}
if corr_plot:
    correlation_matrix = {exp: nf.corrMatrix(raw_data[exp], z_trans=False) for exp in raw_data.keys()}

    threshold_matrix = {exp: nf.significanceCheck(matrix[2], matrix[0], alpha=1, results_folder = results_folder, title = exp,
                                                  names=nodes[exp], plot = corr_plot, include_Negs=False, Anatomy= ROIs,
                                                  savefig = True) \
                        for exp, matrix in correlation_matrix.items()}
correlation_matrix ={exp: nf.corrMatrix(raw_data[exp], z_trans=True) for exp in raw_data.keys()}

threshold_matrix = {exp: nf.significanceCheck(matrix[2], matrix[0],
                                              alpha = 1, names = nodes[exp], include_Negs = False)\
                    for exp, matrix in correlation_matrix.items()}

per = {exp: {str(p): nf.percentile(matrix, p) for p in percentile}\
                    for exp, matrix in threshold_matrix.items()}
#threshold_matrix_an = {exp: nf.significanceCheck(matrix[2], matrix[0],
#                                              alpha = 1, names = nodes[exp], include_Negs = True, Anatomy = ROIs)\
#                    for exp, matrix in correlation_matrix.items()}
graph = {exp: {p: nf.networx(matrix, nodes[exp])[0] for p, matrix in pe.items()} for exp, pe in per.items()}
pos = {exp: {p: nf.networx(matrix, nodes[exp])[1] for p, matrix in pe.items()} for exp, pe in per.items()}


#CLUSTERING
hc_dict = {exp:{p: al.hierarch_clust(matrix, nodes[exp],ROIs.values(),plot = False) for p, matrix in pe.items()} for exp, pe in graph.items()}

for p in percentile:
    p = str(p)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(hc_dict['KET30'][p][0]['Distance Cut'],hc_dict['KET30'][p][0]['Number of Clusters'],c='orange',label='KET30' + p)
    ax1.scatter(hc_dict['SAL'][p][0]['Distance Cut'],hc_dict['SAL'][p][0]['Number of Clusters'],c='grey',label='SAL' + p)
    plt.title(p)
    plt.xlabel("Cut Distance Along Dendrogram")
    plt.ylabel("Number of Clusters")
    plt.legend(loc='upper right')
    plt.show()
    fig.savefig(os.path.join(results_folder, 'HC_Cluster oer dendrogram' +"_"+ p+\
                             '.png'))

color_dict_hc = {exp:{p: plu.grab_color_attributes(matrix[2], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}
pos_dict_hc = {exp:{p: plu.get_position_data(matrix[2], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}

{exp: {p: plu.graph_network(matrix, list(color_dict_hc[exp][p].values()), pos_dict_hc[exp][p], '_'.join(['HC',exp, p]), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}
my_allen_colors = plu.get_allen_colors('ROIs.csv')
{exp: {p: plu.graph_network(matrix, my_allen_colors, pos_dict_hc[exp][p], '_'.join(['HC', exp, p, 'Anatomic']), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}


markov_dict = {exp: {p:al.markov(matrix, nodes[exp]) for p, matrix in pe.items()} for exp, pe in graph.items()}
color_dict = {exp:{p: plu.grab_color_attributes(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in markov_dict.items()}
pos_dict = {exp:{p: plu.get_position_data(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in markov_dict.items()}

{exp: {p: plu.graph_network(matrix, list(color_dict[exp][p].values()), pos_dict[exp][p], '_'.join([exp, p]), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}
my_allen_colors = plu.get_allen_colors('ROIs.csv')
{exp: {p: plu.graph_network(matrix, my_allen_colors, pos_dict[exp][p], '_'.join([exp, p, 'Anatomic']), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}



louvain_dict = {exp: {p:al.louvain(matrix, nodes[e], 100) for p, matrix in pe.items()} for exp, pe in graph.items()}

color_dict_l = {exp:{p: plu.grab_color_attributes(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
pos_dict_l = {exp:{p: plu.get_position_data(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}



{exp: {p: plu.graph_network(matrix, list(color_dict_l[exp][p].values()), pos_dict_l[exp][p], '_'.join(['Louvain',exp, p]), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}
my_allen_colors = plu.get_allen_colors('ROIs.csv')
{exp: {p: plu.graph_network(matrix, my_allen_colors, pos_dict_l[exp][p], '_'.join(['Louvain', exp, p, 'Anatomic']), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}



#GLOBAL EFFICIENCY

global_efficiency ={exp: {p: nx.global_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}
local_efficiency ={exp: {p:nx.local_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}

delta_global_eff = {exp: {p:al.in_silico_deletion(matrix, '_'.join([exp, p]), results_folder, plot =True, save = True) for p, matrix in pe.items()} for exp, pe in graph.items()}

efficiency = {}

for e in global_efficiency.keys():
    efficiency[e] = {'global_efficiency': [[v for p,v in global_efficiency[e].items()], [p for p,v in global_efficiency[e].items()]],
                     'local_efficiency': [[v for p,v in local_efficiency[e].items()], [p for p,v in local_efficiency[e].items()]]}

plu.plot_efficiency_percentile(efficiency, colors, results_folder, save = True)


#check the modularity at different percentiles

modularity_dict= {e: nf.threshold_simulation(m, 0.05, 0.95, 19, algo='markov') for e, m in threshold_matrix.items()}

#FINDING HUBS

node_attributes = {exp:{p:nf.grab_node_attributes(matrix,use_distance=False,compress_to_df=True)for p, matrix in pe.items()} for exp, pe in graph.items()}
WMDZ_PC_attr = {exp:{p:nf.cluster_attributes(matrix,nodes[exp], markov_dict[exp][p][2], make_df = True)for p, matrix in pe.items()} for exp, pe in graph.items()}



Hubs = {exp:{p:nf.findMyHubs(matrix, ROIs_dict, ROIs) for p, matrix in pe.items()} for exp, pe in graph.items()}

[[df[0].to_csv(os.path.join(results_folder, 'Score_All' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Hubs.items()]

[[df[1].to_csv(os.path.join(results_folder, 'Hubs' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Hubs.items()]
