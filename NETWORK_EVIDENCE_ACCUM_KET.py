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
import pingouin as pg




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

def transpose_df(combined_df):
    trans_df = []
    for ind, col in enumerate(combined_df.columns):
        if col not in ['name', 'condition', 'sex', 'task']:
            dict = {'area' : col,'cfos-density' : list(combined_df[col]), 'name' : list(combined_df['name']),
                    'condition': list(combined_df['condition']), 'sex': list(combined_df['sex']),
                    'task': list(combined_df['task'])}
            trans_df.append(pd.DataFrame.from_dict(dict))
    trans_df = pd.concat(trans_df)
    return trans_df
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


# %%
# calc shortest path cause for some reason its having trouble with it inside the function
def shortest(G, sigData):
    short = nx.floyd_warshall_numpy(G, weight='weight')
    shortavg = np.mean(short, axis=0)
    keys = sigData.index.values.tolist()
    vals = shortavg.tolist()
    zip_iterator = zip(keys, vals)
    short_dictionary = dict(zip_iterator)

    return short_dictionary


def FindMyHubs(G, sigData):
    G_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in
                       G.edges(data='weight')}  # creates a dict of calculted distance between all nodes
    nx.set_edge_attributes(G, values=G_distance_dict, name='distance')  # sets the distance as an atribute to all nodes

    cluster_coefficient = nx.clustering(G, weight='weight')  # calculated different measures of importance
    degree_cent = nx.degree_centrality(G)
    eigen_cent = nx.eigenvector_centrality(G, max_iter=100000, weight='weight')
    betweenness = nx.betweenness_centrality(G, weight='distance', normalized=True)
    closeness = nx.closeness_centrality(G, distance='distance')
    degree = list(G.degree(G, weight='weight'))
    communicability = nx.communicability_betweenness_centrality(G)
    short_avg = shortest(G, sigData)

    dict = {'eigen_cent': eigen_cent, 'betweenness': betweenness, 'closeness': closeness,
            'clustering': cluster_coefficient, 'degree_cent': degree_cent, 'communicability': communicability,
            'short_avg': short_avg}

    Results = pd.DataFrame(dict)  # create a frame with all the measures of importance for every region

    # Van den Huevel(2010) - https://www.jneurosci.org/content/30/47/15915
    # used the top or bottom quartiles to determine the hubness of all nodes so here we calculate that.

    eig_quart_upper = Results.eigen_cent.quantile(0.75)  # eig centrality upper and lower quartiles
    # eig_quart_lower = Results.eigen_cent.quantile(.25)
    degree_cent_quart_upper = Results.degree_cent.quantile(0.75)  # weighted degree centrality upper quartiles
    betweeness_quart_upper = Results.betweenness.quantile(0.75)  # betweeness centrality upper quartile
    closeness_quart_upper = Results.closeness.quantile(.75)  # closeness centrality upper quartile
    clustering_quart_lower = Results.clustering.quantile(.25)  # clustering lower quartile
    com_quart_upper = Results.communicability.quantile(.75)  # communicability upper quartile
    short_quart_lower = Results.short_avg.quantile(.25)
    quartiles_dict = {'eing_cent': eig_quart_upper, 'betweeness': betweeness_quart_upper,
                      'clustering': clustering_quart_lower, 'shortness': short_quart_lower}
    # for each significant measure an ROI has add one in the score column, a score >= 2 is considered a hub node.
    Results['score'] = 0
    Results['score'] = np.where((Results['eigen_cent'] > eig_quart_upper), Results['score'] + 1, Results['score'])
    # Results['score'] = np.where((Results['eigen_cent'] < eig_quart_lower), Results['score']+1, Results['score'])
    Results['score'] = np.where((Results['betweenness'] >= betweeness_quart_upper), Results['score'] + 1,
                                Results['score'])
    Results['score'] = np.where((Results['clustering'] <= clustering_quart_lower), Results['score'] + 1,
                                Results['score'])
    # Results['score'] = np.where((Results['communicability'] >= com_quart_upper), Results['score']+1, Results['score'])
    Results['score'] = np.where((Results['short_avg'] <= short_quart_lower), Results['score'] + 1, Results['score'])

    NonHubs = Results[(Results['score'] < 2)].index  # create an index of rois with a score of less than 2 in hubness

    Hubs = Results.drop(NonHubs, errors='ignore')  # create a new frame with only the important nodes

    return Results, Hubs, quartiles_dict

bilateral = False
ROI = False
processed = False
mean = True
percentage_of_change = False
corr_plot = False
persex  =False

percen = np.arange(0.1, 0.7, 0.1)
percent = [round(n, 2) for n in percen]
groups_compare = [['SAL-HC', 'SAL-EAT'], ['SAL-HC', 'KET30-HC'], ['SAL-EAT', 'KET30-EAT'], ['KET30-HC', 'KET30-EAT']]
used_P = '0.25'

extra_name = 'EAT'

filepath = 'C:/Users/cds4/Desktop/Network/EVIDENCE_ACCUMULATION'

mice = ['MABL6-23', 'MABL6-25', 'MABL6-27','FBL6-29', 'FBL6-31','FBL6-33', 'MABL6-36', 'MABL6-38','MABL6-40', 'FBL6-42', 'FBL6-44','FBL6-46', 'MABL6-47', 'MABL6-49', 'FBL6-51', 'MABL6-24', 'MABL6-26', 'MABL6-28','FBL6-30', 'FBL6-32','FBL6-34', 'MABL6-35', 'MABL6-37','MABL6-39', 'FBL6-41', 'FBL6-43','FBL6-45', 'MABL6-48', 'FBL6-50', 'FBL6-52']



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

Allen_areas_dict_filepath = "C:/Users/cds4/Documents/GitHub/network_analysis/Allen_Areas_dict.pickle"
results_folder = 'C:/Users/cds4/Desktop/Network/RESULTS_EAT'
filepath = 'C:/Users/cds4/Desktop/Network/COMBINED'
if ROI:
    data_wrangled = {'KET30-HC': 'C:/Users/cds4/Desktop/Network/COMBINED/KET30_Large_Network_ROI.csv',
                     'SAL-HC': 'C:/Users/cds4/Desktop/Network/COMBINED/SAL_Large_Network_ROI.csv',
                     'KET30-EAT': 'C:/Users/cds4/Desktop/Network/EVIDENCE_ACCUMULATION/KET30_EAT_Large_Network_ROI.csv',
                     'SAL-EAT': 'C:/Users/cds4/Desktop/Network/EVIDENCE_ACCUMULATION/SAL_EAT_Large_Network_ROI.csv'}
else:
    data_wrangled = {'KET30-HC': 'C:/Users/cds4/Desktop/Network/COMBINED/KET30_Large_Network.csv',
                     'SAL-HC': 'C:/Users/cds4/Desktop/Network/COMBINED/SAL_Large_Network.csv',
                     'KET30-EAT': 'C:/Users/cds4/Desktop/Network/EVIDENCE_ACCUMULATION/KET30_EAT_Large_Network.csv',
                     'SAL-EAT': 'C:/Users/cds4/Desktop/Network/EVIDENCE_ACCUMULATION/SAL_EAT_Large_Network.csv'
                     }
animal_dict = {'KET30-HC': ['M11', 'M1', 'M3', 'M5', 'M7', 'M9', 'F11', 'F13', 'F15', 'F1', 'F3', 'F5', 'F7', 'F9'],
               'SAL-HC': ['M10', 'M12', 'M2', 'M4', 'M6', 'M8', 'F10', 'F12', 'F14', 'F2', 'F4', 'F6', 'F8'],
               'KET30-EAT': ['MABL6-40', 'FBL6-42', 'MABL6-47', 'MABL6-25', 'FBL6-33', 'FBL6-46', 'MABL6-49', 'FBL6-51'],
               'SAL-EAT': ['MABL6-48', 'FBL6-50', 'MABL6-39', 'FBL6-43', 'FBL6-34', 'MABL6-28', 'MABL6-35', 'FBL6-52']}

colors_exp = {'KET30-EAT': '#FF6103', 'KET30-HC': '#ED9121', 'SAL-HC': '#030303', 'SAL-EAT': '#A9A9A9'}

color_list = ['#030303', '#A9A9A9', '#ED9121', '#FF6103']



ROIs, Allen_Groups = dw.ag_pickle(Allen_areas_dict_filepath)
ROIs_dict = dw.get_ROIs_dict(filepath)
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
nodes = {exp: nf.loadData(data_wrangled[exp], csv = True)[1] for exp in data_wrangled.keys()}


if percentage_of_change:
    combined_df = combine_exp(raw_data, animal_dict)
    combined_data = combined_df.drop(['condition', 'name'], axis=1)
    combined_df['sex'] = ['Male' if 'M' in n else 'Female' if 'F' in n  else None for n in combined_df.name]
    combined_df['task'] = [exp.split('-')[1] for exp in combined_df.condition]
    #Areas ordered anatomically
    data_mean = {}
    for exp in set(list(combined_df.condition)):
        data_concat = []
        data = combined_df.query('condition==@exp')
        data = data.fillna(0)
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
        trans_df = transpose_df(combined_df)
        trans_df['complete_name'] = [ROIs_dict[area] for area in list(trans_df.area)]
        trans_df['region'] = [ROIs[area] for area in list(trans_df.area)]

        g = sns.catplot(data=trans_df, x="cfos-density", y="area", hue="condition",
                    hue_order = ['SAL-HC', 'SAL-EAT', 'KET30-HC', 'KET30-EAT'], palette = color_list,
                    kind= 'boxen', height =30, aspect=2, orient = 'h')
        g.set_xticklabels(fontsize=30)
        title = 'Change between conditions'

        g.savefig(os.path.join(results_folder, title + \
                                 '.png'), dpi=300)


        region_df = []
        done_animals = []
        for mouse, region in zip(trans_df.name,trans_df.region):
            if '-'.join([mouse, region]) in done_animals:
                continue
            else:
                done_animals.append('-'.join([mouse, region]))
            data = trans_df.query('region == @region and name == @mouse')
            mean = data['cfos-density'].mean()
            condition = list(data['condition'])[0]
            sex = list(data['sex'])[0]
            task = list(data['task'])[0]
            region_df.append([mouse, condition, sex, task, region, mean])
        region_df = pd.DataFrame(region_df, columns=['name', 'condition', 'sex', 'task','region', 'cfos-density'])


        fig, ax  =plt.subplots(nrows = 1, ncols =1,figsize=(12,10))
        sns.barplot(data=region_df, x="region", y="cfos-density", hue="condition", order=['Cerebellum','Cortical plate','Cortical subplate',
                        'Hypothalamus', 'Medulla', 'Midbrain','Pallidum', 'Pons', 'Striatum','Thalamus'],
                        hue_order=['SAL-HC', 'SAL-EAT', 'KET30-HC', 'KET30-EAT'],palette=color_list, ax = ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        title = 'Change between conditions region'

        fig.savefig(os.path.join(results_folder, title + \
                               '.png'), dpi=300)
        stats = []
        posthocs= []
        sig_areas = []
        sig_posthocs = []
        stats_region = []
        posthocs_region = []
        sig_areas_region = []
        sig_posthocs_region = []
        tstud_task = []
        significant_areas = []
        for area in set(trans_df.area):
            data = trans_df.query('area == @area')
            tw = pg.anova(data=data, dv='cfos-density', between=['task', 'condition'])
            tw['area'] = area
            tw['complete_name'] = ROIs_dict[area]
            tw['region'] =ROIs[area]
            stats.append(tw)
            pval_con = list(tw.query('Source == "condition"')['p-unc'])[0]

            ph = pg.pairwise_ttests(dv='cfos-density', between=['task', 'condition'],
                                    data=data, parametric=True, marginal=True,
                                    alpha=0.05, padjust='bonf', effsize='eta-square',
                                    correction='auto',
                                    nan_policy='listwise', return_desc=False, interaction=False,
                                    within_first=False)

            ph['area'] = area
            ph['complete_name'] = ROIs_dict[area]
            ph['region'] = ROIs[area]
            if pval_con < 0.05:
                sig_areas.append(area)
                sig_posthocs.append(ph)

            posthocs.append(ph)
            x = np.array(data.query('task == "EAT"')['cfos-density'].values)
            y = np.array(data.query('task == "HC"')['cfos-density'].values)
            n = np.concatenate([x, y])
            normality = pg.normality(n)
            param = normality['normal'][0]
            if param:
                df = pg.ttest(x, y)
                tstud_task.append([area, df['p-val'][0]])
                if df['p-val'][0] <= 0.05:
                    significant_areas.append(area)
            else:
                df = pg.mwu(x, y)
                tstud_task.append([area, df['p-val'][0]])
                if df['p-val'][0] <= 0.05:
                    significant_areas.append(ROIs_dict[area])

        tstud_task = pd.DataFrame(tstud_task, columns=['area', 'p-value'])
        tstud_task = tstud_task.set_index('area')

        stats = pd.concat(stats)
        posthocs = pd.concat(posthocs)
        sig_posthocs = pd.concat(sig_posthocs)
        stats.to_csv(os.path.join(results_folder, 'Stats_perchange_task' + \
                                  '.csv'))

        posthocs.to_csv(os.path.join(results_folder, 'pothocs_perchange_task' + \
                                     '.csv'))
        sig_posthocs.to_csv(os.path.join(results_folder, 'significant_pothocs_perchange_task' + \
                                     '.csv'))

        for area in set(region_df.region):
            data = region_df.query('region == @area')
            tw = pg.anova(data=data, dv='cfos-density', between=['task', 'condition'])
            tw['region'] = area
            stats_region.append(tw)
            pval_con = list(tw.query('Source == "condition"')['p-unc'])[0]

            ph = pg.pairwise_ttests(dv='cfos-density', between=['task', 'condition'],
                                    data=data, parametric=True, marginal=True,
                                    alpha=0.05, padjust='bonf', effsize='eta-square',
                                    correction='auto',
                                    nan_policy='listwise', return_desc=False, interaction=False,
                                    within_first=False)

            ph['region'] = area

            if pval_con < 0.05:
                sig_areas_region.append(area)
                sig_posthocs_region.append(ph)

            posthocs_region.append(ph)

        stats_region = pd.concat(stats_region)
        posthocs_region = pd.concat(posthocs_region)
        sig_posthocs_region = pd.concat(sig_posthocs_region)
        stats.to_csv(os.path.join(results_folder, 'Stats_region' + \
                                  '.csv'))

        posthocs.to_csv(os.path.join(results_folder, 'pothocs_region' + \
                                     '.csv'))
        sig_posthocs.to_csv(os.path.join(results_folder, 'significant_pothocs_region' + \
                                         '.csv'))

""""
NETWORK ANALYSES
"""

percen = np.arange(0.25, 1.1, 0.25)
percent = [round(n, 2) for n in percen]
raw_data = {exp: nf.loadData(data_wrangled[exp], csv = True)[0] for exp in data_wrangled.keys()}
correlation_matrix ={exp: nf.corrMatrix(raw_data[exp], z_trans=True) for exp in raw_data.keys()}


threshold_matrix = {exp: significanceCheck(matrix[2], matrix[0],
                                              alpha = 1, names = nodes[exp], title = exp, results_folder=results_folder, include_Negs = False)\
                    for exp, matrix in correlation_matrix.items()}

per_adj = {exp: {str(p): percentile(matrix[0], p, names = nodes[exp], title = '_'.join([exp, str(p)]), results_folder = results_folder, plot = True,Anatomy= ROIs) for p in percent}\
                    for exp, matrix in threshold_matrix.items()}


graph = {exp: {p: nf.networx(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}


#GLOBAL EFFICIENCY

global_efficiency ={exp: {p: nx.global_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}
local_efficiency ={exp: {p:nx.local_efficiency(matrix) for p, matrix in pe.items()} for exp, pe in graph.items()}

#delta_global_eff = {exp: {p:al.in_silico_deletion(matrix, '_'.join([exp, p]), results_folder, plot =True, save = True) for p, matrix in pe.items()} for exp, pe in graph.items()}

efficiency = {}

for e in global_efficiency.keys():
    efficiency[e] = {'global_efficiency': [[v for p,v in global_efficiency[e].items()], [p for p,v in global_efficiency[e].items()]],
                     'local_efficiency': [[v for p,v in local_efficiency[e].items()], [p for p,v in local_efficiency[e].items()]]}

plu.plot_efficiency_percentile(efficiency, colors_exp, results_folder, save = True)

### DRAW NETWORKS

louvain_dict = {exp: {p:al.louvain(matrix, nodes[e], 100) for p, matrix in pe.items()} for exp, pe in graph.items()}

color_dict_l = {exp:{p: plu.grab_color_attributes(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}
pos_dict_l = {exp:{p: plu.get_position_data(matrix[0], nodes[exp]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}



{exp: {p: plu.graph_network(matrix, list(color_dict_l[exp][p].values()), pos_dict_l[exp][p], '_'.join(['Louvain',exp, p]), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}
allen_list_alphabetical = sorted(colors.keys())
allen_colors = [color for color in sns.color_palette('Paired', len(allen_list_alphabetical))]
allen_color_dict = {group: color for group, color in zip(allen_list_alphabetical, allen_colors)}
color_list_allen = []
for area in ROIs:
    color_list_allen.append(allen_color_dict[ROIs[area]])

{exp: {p: plu.graph_network(matrix, color_list_allen, pos_dict_l[exp][p], '_'.join(['Louvain', exp, p, 'Anatomic']), results_folder, save = True)\
       for p, matrix in pe.items()} for exp, pe in graph.items()}



results = {exp: {p: FindMyHubs(matrix , per_adj[exp][p][1]) for p, matrix in pe.items()} for exp, pe in graph.items()}
dw.pickle_save(results, results_folder + '/results.pickle')
Modules= {exp: {p: np.array([str(x) for x in matrix[3].flatten()]) for p, matrix in pe.items()} for exp, pe in louvain_dict.items()}

index = {exp: {p: matrix[1].index.copy() for p, matrix in pe.items()} for exp, pe in results.items()}

indexx_full_network = {exp: {p: matrix[0].index.copy() for p, matrix in pe.items()} for exp, pe in results.items()}

modules_full_network = {exp: {p: pd.DataFrame(matrix, columns=['Module'],index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in Modules.items()}

within_module_degree_FN = {exp: {p:  brainconn.centrality.module_degree_zscore(matrix[0], modules_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}

FN = {exp: {p:  brainconn.centrality.participation_coef_sign(matrix[0], modules_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in per_adj.items()}


Z_data_FN = {exp: {p:  pd.DataFrame(matrix, columns=['Within module degree'],index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in within_module_degree_FN.items()}
P_data_FN = {exp: {p:  pd.DataFrame(matrix[0], columns=['Participation coefficient'], index = indexx_full_network[exp][p]) for p, matrix in pe.items()} for exp, pe in FN.items()}



Network_all_data = {exp: {p: pd.concat([matrix[0], Z_data_FN[exp][p], P_data_FN[exp][p], modules_full_network[exp][p]], axis=1) for p, matrix in pe.items()} for exp, pe in results.items()}



Hubs_per_area = {}
Network_area_combined = {}
for exp, pe in Network_all_data.items():
    for p,m in pe.items():
        m['Complete name'] = [ROIs_dict[getattr(r, 'Index')] for r in m.itertuples()]
        m['Region'] = [ROIs[getattr(r, 'Index')] for r in m.itertuples()]
        Network_all_data[exp][p] = m[['Complete name', 'Region', 'score', 'Module',  'Within module degree',
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

data_stats['task'] = [x.split('-')[1] for x in data_stats.exp]
stats_network = []
posthocs_network = []
sig_variables = []
sig_posthocs_network = []
for col in list(data_stats.columns):
    if col not in no_stats:
        n = np.array(data_stats[col])
        normality = pg.normality(n)
        param = normality['normal'][0]
        tw = pg.anova(data=data_stats, dv=col, between=['task', 'exp'])
        tw['variable'] = col
        stats_network.append(tw)

        # pval_con = list(tw.query('Source == "exp"')['p-unc'])[0]

        ph = pg.pairwise_ttests(dv=col, between=['task', 'exp'],
                                data=data_stats, parametric=True, marginal=True,
                                alpha=0.05, padjust='bonf', effsize='eta-square',
                                correction='auto',
                                nan_policy='listwise', return_desc=False, interaction=False,
                                within_first=False)

        ph['variable'] = col
        # if pval_con < 0.05:
        #     sig_variables.append(col)
        #     sig_posthocs_network.append(ph)

        posthocs_network.append(ph)


stats_network = pd.concat(stats_network)
posthocs_network = pd.concat(posthocs_network)
sig_posthocs_network = pd.concat(sig_posthocs_network)

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


fig.legend(handles=plu.create_custom_legend(colors_exp), prop={"size": 25}, loc = 4)
# plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
fig.savefig(os.path.join(results_folder, 'Number of hubs per area' + exp + '_' + used_P + '.png'), dpi=300)

[[df.to_csv(os.path.join(results_folder, 'Score_All' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Network_all_data.items()]

[[df.rename(columns = {'Module': 'Modularity'}).to_csv(os.path.join(results_folder, 'gephi' + e +"_"+ p+\
                             '.csv')) for p, df in pe.items()] for e, pe in Network_all_data.items()]


#Need to be done two compare across conditions
shared_hubs = {}
hubs_exp = {}
different_hubs = {}
for exp,per in Network_all_data.items():
    hubs_exp[exp] = list(per[used_P][per[used_P].score > 1].index)

for group in groups_compare:
    shared_hubs['_'.join(group)] = set(hubs_exp[group[0]]).intersection(hubs_exp[group[1]])
    different_hubs['_'.join(group)] = {group[0] : [ element for element in hubs_exp[group[0]] if element not in hubs_exp[group[1]]],\
                                      group[1]: [element for element in hubs_exp[group[1]] if element not in hubs_exp[group[0]]]}

"""
UMAP
"""

colors_list =[len(animal_dict[exp]) * [colors_exp[exp]] for exp in colors_exp.keys()]
color_list = colors_list[0] + colors_list[1] + colors_list[2] + colors_list[3]



combined_df = combine_exp(raw_data, animal_dict)
combined_data = combined_df.drop(['condition', 'name'], axis = 1)

n_neighbors_list = np.arange(2,14,1)

min_dist_list =(0.0, 0.1, 0.25, 0.5, 0.8, 0.99)

components_list =[1,2,3]

components = 1
neighbors = 7
min_dist = 0.5
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




