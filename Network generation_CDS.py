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
<<<<<<< Updated upstream:Network generation_CDS.py
bilateral = False
ROI = False
=======
import itertools
import networkx as nx
import pandas as pd
import numpy as np
import scipy.special as sc
import seaborn as sns
import matplotlib.patches as mpatches
from bokeh.plotting import from_networkx
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.io import output_file, show
from bokeh.models.graphs import from_networkx
from bokeh.models import (BoxZoomTool, Circle, HoverTool,
                          MultiLine, Plot, Range1d, ResetTool,)
from bokeh.palettes import Spectral4
from matplotlib import pyplot as plt
import brainconn
import Data_wrangling as dw
import network_analysis.NetworkFunctions as nf
import network_analysis.plotting_utils as plu
import network_analysis.algorithms as al
import os
import umap
from sklearn.preprocessing import StandardScaler
import seaborn as sns

def combine_exp(raw_data,animal_dict):
    combined_df = []
    for e,df in raw_data.items():
        df['condition'] = e
        df['name'] = animal_dict[e]
        combined_df.append(df)
    combined_df = pd.concat(combined_df)
    return combined_df

def draw_umap(data, mice_list, n_neighbors=5, min_dist=0.01, n_components=2, metric='euclidean', title='', results_folder = '', cmap = 'Spectral', save = False,):
    fit = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric
    )
    u = fit.fit_transform(data);
    fig = plt.figure()
    if cmap == 'Spectral':
        if n_components == 1:
            ax = fig.add_subplot(111)
            sp = ax.scatter(u[:,0], range(len(u)), c=[x for x, n in enumerate(mice_list)], cmap=cmap)
            plt.gca().set_aspect('equal', 'datalim')
            fig.colorbar(sp, boundaries=np.arange(len(mice_list)+1) - 0.5).set_ticks(ticks=\
                                                                 np.arange(len(mice_list)),
                                                                 labels=mice_list)
        if n_components == 2:
            ax = fig.add_subplot(111)
            sp = ax.scatter(u[:,0], u[:,1], c=[x for x, n in enumerate(mice_list)], cmap=cmap)
            plt.gca().set_aspect('equal', 'datalim')
            fig.colorbar(sp, boundaries=np.arange(len(mice_list)+1) - 0.5).set_ticks(ticks=\
                                                                 np.arange(len(mice_list)),
                                                                 labels=mice_list)
        if n_components == 3:
            ax = fig.add_subplot(111, projection='3d')
            sp = ax.scatter(u[:,0], u[:,1], u[:,2], c=[x for x, n in enumerate(mice_list)], cmap=cmap,
                       s=100)
            plt.gca().set_aspect('auto', 'datalim')
            plt.colorbar(sp, boundaries=np.arange(len(mice_list)+1) - 0.5).set_ticks(ticks=\
                                                                 np.arange(len(mice_list)),
                                                                 labels=mice_list)
    else:
        if n_components == 1:
            ax = fig.add_subplot(111)
            sp = ax.scatter(u[:,0], range(len(u)), c=cmap)
            plt.gca().set_aspect('equal', 'datalim')
        if n_components == 2:
            ax = fig.add_subplot(111)
            sp = ax.scatter(u[:,0], u[:,1], c=cmap)
            plt.gca().set_aspect('equal', 'datalim')

        if n_components == 3:
            ax = fig.add_subplot(111, projection='3d')
            sp = ax.scatter(u[:,0], u[:,1], u[:,2], c=cmap,
                       s=100)
            plt.gca().set_aspect('auto', 'datalim')
    plt.title(title, fontsize=18)
    if save:
        fig.savefig(os.path.join(results_folder,title+\
                                 '.png'))
    return u



bilateral = False
ROI = True
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
processed = False
mean = True
percentage_of_change = False
corr_plot = False
<<<<<<< Updated upstream:Network generation_CDS.py
percentile = np.arange(0.05, 1, 0.05)
percentile = [round(n, 2) for n in percentile]
#percentile = [0.2]

colors = {'KET30': 'orange', 'SAL': 'grey'}
=======
persex  =False
percen = np.arange(0.1, 0.7, 0.1)
percent = [round(n, 2) for n in percen]
#percentile = [0.2]
sex = 'C'

colors_exp = {'KET30': 'orange', 'SAL': 'grey'}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

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

<<<<<<< Updated upstream:Network generation_CDS.py

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
=======
Allen_areas_dict_filepath = "C:/Users/cds4/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle"


if sex == 'F':
    filepath = 'C:/Users/cds4/Desktop/Network/NEW_DATA'
    results_folder = 'C:/Users/cds4/Desktop/Network/RESULTS_FEM'
    data_wrangled, filepath, results_folder, data_mean = dw.dataManager(filepath, results_folder, sex, bilateral, ROI, processed, mean, new = True)
    if ROI:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/NEW_DATA/KET30_Large_Network_ROI.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/NEW_DATA/SAL_Large_Network_ROI.csv'}
    else:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/NEW_DATA/KET30_Large_Network.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/NEW_DATA/SAL_Large_Network.csv'}
elif sex == 'M':
    filepath = 'C:/Users/cds4/Desktop/Network/MALE'
    results_folder = 'C:/Users/cds4/Desktop/Network/Results'
    data_wrangled, filepath, results_folder, data_mean = dw.dataManager(filepath, results_folder, sex,bilateral, ROI, processed, mean)
    if ROI:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/MALE/KET30_Large_Network_ROI.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/MALE/SAL_Large_Network_ROI.csv'}
    else:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/MALE/KET30_Large_Network.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/MALE/SAL_Large_Network.csv'}
else:
    filepath = 'C:/Users/cds4/Desktop/Network/COMBINED'
    results_folder = 'C:/Users/cds4/Desktop/Network/RESULTS_COMBINED'
    #data_wrangled, filepath, results_folder, data_mean = dw.dataManager(filepath, results_folder, sex, bilateral, ROI,
    #                                                                    processed, mean)
    if ROI:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/COMBINED/KET30_Large_Network_ROI.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/COMBINED/SAL_Large_Network_ROI.csv'}
    else:
        data_wrangled = {'KET30': 'C:/Users/cds4/Desktop/Network/COMBINED/KET30_Large_Network.csv',
                         'SAL': 'C:/Users/cds4/Desktop/Network/COMBINED/SAL_Large_Network.csv'}
    persex = True
if sex == 'F':
    animal_dict = {'KET30': ['F11', 'F13', 'F15', 'F1', 'F3', 'F5', 'F7', 'F9', ],
                   'SAL': ['F10', 'F12', 'F14', 'F2', 'F4', 'F6', 'F8']}
elif sex == 'M':
    animal_dict = {'KET30': ['M11', 'M1', 'M3', 'M5', 'M7', 'M9', ],
                   'SAL': ['M10', 'M12', 'M2', 'M4', 'M6', 'M8',]}
else:
    animal_dict = {'KET30': ['M11','M1', 'M3', 'M5', 'M7', 'M9', 'F11', 'F13', 'F15','F1', 'F3', 'F5', 'F7', 'F9',],
              'SAL': ['M10', 'M12','M2', 'M4', 'M6', 'M8', 'F10', 'F12', 'F14','F2', 'F4', 'F6', 'F8']}
    animal_dict1 = {'KET30_F': ['F11', 'F13', 'F15', 'F1', 'F3', 'F5', 'F7', 'F9', ],
                   'SAL_F': ['F10', 'F12', 'F14', 'F2', 'F4', 'F6', 'F8'],
                   'KET30_M': ['M11', 'M1', 'M3', 'M5', 'M7', 'M9', ],
                   'SAL_M': ['M10', 'M12', 'M2', 'M4', 'M6', 'M8', ]
                   }
    colors_exp = {'KET30_M': 'orange', 'KET30_F':'#ED9121', 'SAL_M': 'darkgrey','SAL_F': 'grey'}


#extra_discard = ['MB', 'MY', 'P']
#for e in data_wrangled.keys():
 #   data_wrangled[e] = data_wrangled[e].drop(extra_discard, axis = 1)
# if ROI:
#     data_wrangled = {'KET30':  'C:/Users/cds4/Desktop/Network/KET30_Large_Network_ROI.csv',
#                  'SAL':'C:/Users/cds4/Desktop/Network/SAL_Large_Network_ROI.csv' }
# else:
#     data_wrangled = {'KET30':  'C:/Users/cds4/Desktop/Network/KET30_Large_Network.csv',
#                  'SAL':'C:/Users/cds4/Desktop/Network/SAL_Large_Network.csv' }

ROIs, Allen_Groups = dw.ag_pickle(Allen_areas_dict_filepath)
ROIs_dict = dw.get_ROIs_dict(filepath)
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
nodes = {exp: nf.loadData(data_wrangled[exp], csv = True)[1] for exp in data_wrangled.keys()}


if percentage_of_change:
    combined_df = combine_exp(raw_data, animal_dict)
    combined_data = combined_df.drop(['condition', 'name'], axis=1)
    combined_df['sex'] = ['Male' if 'M' in n else 'Female' if 'F' in n  else None for n in combined_df.name]
    #Areas ordered anatomically
    if sex == 'C':
        data_mean = {}
        for exp in set(list(combined_df.condition)):
            data_concat = []
            data = combined_df.query('condition==@exp')
            data= data.fillna(0)
            dm = pd.DataFrame([data.mean().index, data.mean().values])
            cols = dm.iloc[0]
            dm = dm[1:]
            dm.columns = cols
            if persex:
                dm['comparison'] = 'both'
                data_concat.append(dm)
                for s in set(list(data.sex)):
                    data_sex = data.query('sex==@s')
                    df = pd.DataFrame([data_sex.mean().index, data_sex.mean().values])
                    cols = df.iloc[0]
                    df = df[1:]
                    df.columns = cols
                    df['comparison'] = s
                    data_concat.append(df)
                data_mean[exp] = pd.concat(data_concat)
            else:
                data_mean[exp] = dm

    norm_group = 'SAL'
    if persex:
        order_areas = data_mean['SAL'].drop('comparison', axis =1).columns

    else:
        order_areas = data_mean['SAL'].columns
    order_complete_name = [ROIs_dict[ab] for ab in order_areas]
    order_br = [ROIs[a] for a in order_areas]
    colors = plu.color_allens(sorted(set(order_br)))
    color_order = [colors[b] for b in order_br]
    #percentage of change with respect of the baseline
    percentage_change_dict = q.percentage_change(data_mean, norm_group = 'SAL', persex = persex)
    if persex:
        stats_dict_anova =  q.twanova(combined_df, threshold = 0.05)
        stats_dict = {e: q.tstud(raw_data['KET30'], raw_data[norm_group], threshold = 0.05) for e, data in raw_data.items() if e != norm_group}
    else:
        stats_dict = {e: q.tstud(raw_data['KET30'], raw_data[norm_group], threshold = 0.05) for e, data in raw_data.items() if e != norm_group}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
    pc_df = {}
    sig_pc_df = {}
    sig_pc_dict = {}
    order_sig_br = {}
<<<<<<< Updated upstream:Network generation_CDS.py
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
=======
    for e in stats_dict.keys():
        if sex =='C':
            exp = '_'.join([e, 'both'])
        else:
            exp = e
        pc_df[e] = pd.concat([percentage_change_dict[exp][2], stats_dict[e][0]], axis=1)
        cn = []
        reg = []
        sig_df = []
        acr = []
        for r in pc_df[e].iterrows():
            cn.append(ROIs_dict[r[0]])
            reg.append(ROIs[r[0]])
            acr.append(r[0])
        pc_df[e]['Complete name'] = cn
        pc_df[e]['Region'] = reg
        pc_df[e]['Acronym'] = acr
        pc_df[e] = pc_df[e][['Acronym', 'Complete name', 'Region','Percentage of change', 'p-value']]
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
        pc_df[e] = pc_df[e].sort_values('Region')
        sig_df = pc_df[e].loc[stats_dict[e][1]]
        sig_values = np.array(sig_df['Percentage of change'])
        reg_str = np.array(sig_df['Region'])
<<<<<<< Updated upstream:Network generation_CDS.py
        colors_sig = [colors[r] for r in list(sig_df['Region'])]
        order_sig_br[e] = [np.arange(len(list(sig_df['Region']))),list(sig_df['Complete name'])]

        sig_pc_df[e] = [sig_values, reg_str, colors_sig, sig_df]
=======
        colors_sig = {r:colors[r] for r in list(sig_df['Region'])}
        colors_sig2 = [colors[r] for r in list(sig_df['Region'])]
        order_sig_br[e] = [np.arange(len(list(sig_df['Region']))),list(sig_df['Complete name']), list(sig_df['Acronym'])]

        sig_pc_df[e] = [sig_values, reg_str, colors_sig, sig_df,colors_sig2]
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

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


<<<<<<< Updated upstream:Network generation_CDS.py
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
=======
    fig,ax = plt.subplots(figsize=(34,15))
    for e in sig_pc_df.keys():
        ax.bar(order_sig_br[e][0], sig_pc_df[e][0], width = 0.6, color= sig_pc_df[e][4], tick_label = sig_pc_df[e][1])
    #ax.set_xticklabels(percentage_change_dict[e][1], rotation = 90, fontsize= 30)
    ax.set_xticklabels(order_sig_br[e][2], rotation = 90, fontsize= 20)
    ylabels = ax.get_yticks().tolist()
    label_format = '{:,.0f}'
    ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
    ax.set_yticklabels([label_format.format(x) for x in ylabels],  fontsize= 20)
    fig.legend(handles = plu.create_custom_legend(sig_pc_df[e][2]),  fontsize = 'large', prop={"size":20})
    title = 'Percentage of significant change ' + e
    #ax.set_title(title, fontsize = 70)
    ax.set_ylabel('Percentage of change', fontsize = 25)
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
    fig.savefig(os.path.join(results_folder, title+'.png'))



#Network generation
<<<<<<< Updated upstream:Network generation_CDS.py
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
nodes = {exp: nf.loadData(data_wrangled[exp], csv = True)[1] for exp in data_wrangled.keys()}
=======
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
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

<<<<<<< Updated upstream:Network generation_CDS.py
per = {exp: {str(p): nf.percentile(matrix, p) for p in percentile}\
=======
per_adj = {exp: {str(p): nf.percentile(matrix, p) for p in percent}\
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
                    for exp, matrix in threshold_matrix.items()}
#threshold_matrix_an = {exp: nf.significanceCheck(matrix[2], matrix[0],
#                                              alpha = 1, names = nodes[exp], include_Negs = True, Anatomy = ROIs)\
#                    for exp, matrix in correlation_matrix.items()}
<<<<<<< Updated upstream:Network generation_CDS.py
graph = {exp: {p: nf.networx(matrix, nodes[exp])[0] for p, matrix in pe.items()} for exp, pe in per.items()}
pos = {exp: {p: nf.networx(matrix, nodes[exp])[1] for p, matrix in pe.items()} for exp, pe in per.items()}
=======
graph = {exp: {p: nf.networx(matrix, nodes[exp]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}
#pos = {exp: {p: nf.networx(matrix, nodes[exp])[1] for p, matrix in pe.items()} for exp, pe in per.items()}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py


#CLUSTERING
hc_dict = {exp:{p: al.hierarch_clust(matrix, nodes[exp],ROIs.values(),plot = False) for p, matrix in pe.items()} for exp, pe in graph.items()}

<<<<<<< Updated upstream:Network generation_CDS.py
for p in percentile:
=======
for p in percent:
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
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

<<<<<<< Updated upstream:Network generation_CDS.py
color_dict_hc = {exp:{p: plu.grab_color_attributes(matrix[2], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}
pos_dict_hc = {exp:{p: plu.get_position_data(matrix[2], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}
=======
color_dict_hc = {exp:{p: plu.grab_color_attributes(matrix[3], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}
pos_dict_hc = {exp:{p: plu.get_position_data(matrix[3], nodes[exp]) for p, matrix in pe.items()} for exp, pe in hc_dict.items()}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

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

<<<<<<< Updated upstream:Network generation_CDS.py
color_dict_l = {exp:{p: plu.grab_color_attributes(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
pos_dict_l = {exp:{p: plu.get_position_data(matrix[1], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
=======
color_dict_l = {exp:{p: plu.grab_color_attributes(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
pos_dict_l = {exp:{p: plu.get_position_data(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py



{exp: {p: plu.graph_network(matrix, list(color_dict_l[exp][p].values()), pos_dict_l[exp][p], '_'.join(['Louvain',exp, p]), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}
my_allen_colors = plu.get_allen_colors('ROIs.csv')
{exp: {p: plu.graph_network(matrix, my_allen_colors, pos_dict_l[exp][p], '_'.join(['Louvain', exp, p, 'Anatomic']), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}



#GLOBAL EFFICIENCY

global_efficiency ={exp: {p: nx.global_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}
local_efficiency ={exp: {p:nx.local_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}

<<<<<<< Updated upstream:Network generation_CDS.py
delta_global_eff = {exp: {p:al.in_silico_deletion(matrix, '_'.join([exp, p]), results_folder, plot =True, save = True) for p, matrix in pe.items()} for exp, pe in graph.items()}
=======
#delta_global_eff = {exp: {p:al.in_silico_deletion(matrix, '_'.join([exp, p]), results_folder, plot =True, save = True) for p, matrix in pe.items()} for exp, pe in graph.items()}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

efficiency = {}

for e in global_efficiency.keys():
    efficiency[e] = {'global_efficiency': [[v for p,v in global_efficiency[e].items()], [p for p,v in global_efficiency[e].items()]],
                     'local_efficiency': [[v for p,v in local_efficiency[e].items()], [p for p,v in local_efficiency[e].items()]]}

<<<<<<< Updated upstream:Network generation_CDS.py
plu.plot_efficiency_percentile(efficiency, colors, results_folder, save = True)
=======
plu.plot_efficiency_percentile(efficiency, colors_exp, results_folder, save = True)
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py


#check the modularity at different percentiles

<<<<<<< Updated upstream:Network generation_CDS.py
modularity_dict= {e: nf.threshold_simulation(m, 0.05, 0.95, 19, algo='markov') for e, m in threshold_matrix.items()}
=======
modularity_dict= {e: nf.threshold_simulation(m, nodes[e], 0.1, 0.7, 7, algo='markov') for e, m in threshold_matrix.items()}

plu.plot_efficiency_percentile(modularity_dict, colors, results_folder, save = True)

>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

#FINDING HUBS

node_attributes = {exp:{p:nf.grab_node_attributes(matrix,use_distance=False,compress_to_df=True)for p, matrix in pe.items()} for exp, pe in graph.items()}
WMDZ_PC_attr = {exp:{p:nf.cluster_attributes(matrix,nodes[exp], markov_dict[exp][p][2], make_df = True)for p, matrix in pe.items()} for exp, pe in graph.items()}

<<<<<<< Updated upstream:Network generation_CDS.py


Hubs = {exp:{p:nf.findMyHubs(matrix, ROIs_dict, ROIs) for p, matrix in pe.items()} for exp, pe in graph.items()}
=======
pc_df = pd.read_csv("C:/Users/cds4/Desktop/Network/Results/pc_dfKET30.csv")
pc_df = pc_df.set_index('Unnamed: 0')

Hubs = {exp:{p:nf.findMyHubs(matrix) for p, matrix in pe.items()} for exp, pe in node_attributes.items()}

#Hubs = {exp:{p:nf.findMyHubs(matrix, ROIs_dict, ROIs) for p, matrix in pe.items()} for exp, pe in graph.items()}

Atrr_df = {exp:{p:nf.combine_node_attrs(matrix[0], WMDZ_PC_attr[exp][p], pc_df, global_efficiency[exp][p])for p, matrix in pe.items()} for exp, pe in Hubs.items()}

Hubs_atrr_df = {exp:{p:df.query('Hub_Score >= 3')for p, df in pe.items()} for exp, pe in Atrr_df.items()}

Main_Hubs = {}
hubs_list = {}
for e,per in Hubs_atrr_df.items():
    hubs_list[e] = []
    for p,df in per.items():
        hubs_list[e].append(list(df.index))

Hubs_chain = {e: list(itertools.chain.from_iterable(h)) for e,h in hubs_list.items()}



for e,h in Hubs_chain.items():
    Main_Hubs[e] = {}
    for a in h:
        Main_Hubs[e][a] = h.count(a)

Main_Hubs_df = {e: pd.DataFrame.from_dict(h, orient = 'index').sort_values(0, ascending = False).rename(columns={0 : 'Score_Hubs'}) for e,h in Main_Hubs.items()}

no_hubs_list = {e :[h for h in list(n.values()) if not h in list(Main_Hubs_df[e].index)] for e,n in nodes.items()}

hubs_pc_df = {e :pc_df.drop(nh) for e, nh in no_hubs_list.items()}

Main_Hubs_combined = {e : pd.merge(Main_Hubs_df[e].sort_index(), hubs_pc_df[e].sort_index(), left_index =True, right_index =True)[['Complete name', 'Region', 'Score_Hubs', 'Percentage of change', 'p-value']].sort_values('Score_Hubs', ascending = False) for e in Main_Hubs_df.keys()}
>>>>>>> Stashed changes:development_mats/Network generation_CDS.py

[[df[0].to_csv(os.path.join(results_folder, 'Score_All' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Hubs.items()]

[[df[1].to_csv(os.path.join(results_folder, 'Hubs' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Hubs.items()]
<<<<<<< Updated upstream:Network generation_CDS.py
=======

[[df.to_csv(os.path.join(results_folder, 'Atributes' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Atrr_df.items()]

[df.to_csv(os.path.join(results_folder, 'Score_percentile' + e +\
                             '.csv'))  for e, df in Main_Hubs_combined.items()]




similarity = [al.Similarity(per_adj['KET30'][str(p)], per_adj['SAL'][str(p)]) for p in percent]
dPtargets = {}

for target in Hubs_chain['KET30']:
    dPtargets[target] = {}
    for p in percent:
        dPtargets[target][str(p)] = al.disruptPropagate(graph['KET30'][str(p)], target)

graph_dp = {area: {p: nf.networx(matrix, nodes['KET30']) for p, matrix in pe.items()} for area, pe in dPtargets.items()}


"""
UMAP
"""

colors_list =[len(animal_dict1[exp]) * [colors_exp[exp]] for exp in colors_exp.keys()]
color_list = colors_list[0] + colors_list[1] + colors_list[2] + colors_list[3]



combined_df = combine_exp(raw_data, animal_dict)
combined_data = combined_df.drop(['condition', 'name'], axis = 1)

n_neighbors_list = np.arange(2,14,1)

min_dist_list =(0.0, 0.1, 0.25, 0.5, 0.8, 0.99)

components_list =[1,2,3]

components = 1
neighbors = 4
min_dist = 0.01
metric = 'euclidean'
title = 'UMAP projection of the c-fos activity'
#results_folder = 'C:/Users/cds4/Desktop/Network/Results'


plt.ion()

embedding = draw_umap(combined_data, list(combined_df.name), neighbors, min_dist, components, metric,
          title, results_folder, color_list, save = True)

m_embedding = {str(m): draw_umap(combined_data, list(combined_df.name), neighbors, m, components, metric,
                          title + '_min_dist_' + str(m), results_folder, color_list, save=True) for m in min_dist_list}

n_embedding = {str(n): draw_umap(combined_data, list(combined_df.name), n, min_dist, components, metric,
                          title + '_neighbors_' + str(n), results_folder, color_list, save=True) for n in n_neighbors_list}

c_embedding = {str(c): draw_umap(combined_data, list(combined_df.name), neighbors, min_dist, c, metric,
                          title + '_components_' + str(c), results_folder, color_list, save=True) for c in components_list}


>>>>>>> Stashed changes:development_mats/Network generation_CDS.py
