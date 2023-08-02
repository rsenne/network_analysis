import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os


def grab_color_attributes(cluster_list, node_dict):
    color_list = [color for color in sns.color_palette('colorblind', len(cluster_list))]
    color_dict = {}
    for i in range(0, len(cluster_list)):
        for j in cluster_list[i]:
            color_dict[node_dict[j]] = color_list[i]
    color_dict_sorted = {area: color for area, color in sorted(color_dict.items(), key=lambda ele: ele[0])}
    return color_dict_sorted

def get_allen_colors(allen_rois):
    allen_df = pd.read_csv(allen_rois)
<<<<<<< Updated upstream
    allen_list = list(set(allen_df['Allen Group Name']))
    allen_colors = [color for color in sns.color_palette('colorblind', len(allen_list))]
    allen_color_dict = {group:color for group, color in zip(allen_list, allen_colors)}
    allen_dict = {abbrev:name for abbrev,name in zip(allen_df['Abbreviation'], allen_df['Allen Group Name'])}
=======
    allen_list = list(set(allen_df['Allen Area']))
    allen_list_alphabetical = sorted(allen_list)
    allen_colors = [color for color in sns.color_palette('Paired', len(allen_list_alphabetical))]
    allen_color_dict = {group:color for group, color in zip(allen_list_alphabetical, allen_colors)}
    allen_dict = {abbrev:name for abbrev,name in zip(allen_df['Abbreviation'], allen_df['Allen Area'])}
>>>>>>> Stashed changes
    color_list = []
    for area in allen_dict:
        color_list.append(allen_color_dict[allen_dict[area]])
    return color_list


def sunflower_theta(n):
    golden_ratio = ((1 + 5 ** 0.5) / 2) ** 2
    return 2 * np.pi / golden_ratio * n


def sunflower_r(n, c=1.5):
    return c * (n ** 0.5)


def get_point_cloud(k=0):
    n = [i for i in range(1, k + 1)]
    r = [sunflower_r(i) for i in n]
    theta = [sunflower_theta(j) for j in n]
    point_cloud_x = (r * np.cos(theta))
    point_cloud_y = (r * np.sin(theta))
    point_cloud = [[x, y] for x, y in zip(point_cloud_x, point_cloud_y)]
    return point_cloud


def get_position_data(cluster_list, node_names):
    number_of_clusters = len(cluster_list)
    nodes_list = [x for x in range(0, number_of_clusters)]
    pos_graph = nx.Graph()
    pos_graph.add_nodes_from(nodes_list)
    pos = nx.circular_layout(pos_graph, scale=40, dim=2)
    num_of_nodes = [len(node) for node in cluster_list]
    point_clouds = [get_point_cloud(lens) for lens in num_of_nodes]
    for i in range(len(point_clouds)):
        for j in range(len(point_clouds[i])):
            point_clouds[i][j][0] += pos[i][0]
            point_clouds[i][j][1] += pos[i][1]
    point_cloud_map = {cluster: pos_list for cluster, pos_list in enumerate(point_clouds)}
    pos_dict = {}
    for i in range(len(cluster_list)):
        for j in range(len(cluster_list[i])):
            pos_dict.update({node_names[cluster_list[i][j]]: np.array(point_cloud_map[i][j])})
    nx.rescale_layout_dict(pos_dict, 10)
    return pos_dict


<<<<<<< Updated upstream
def graph_network(G, color_list, pos_dict, title, filepath = '', save = False):
=======
def graph_network(G, color_list, pos_dict, title ='', results_folder ='', save = False):
>>>>>>> Stashed changes
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos_dict, width=1, edge_color=edge_colors, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(G, pos=pos_dict, node_size=node_sizes, node_color=color_list, linewidths=1,
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos=pos_dict)
    ax.margins(0.1, 0.05)
    plt.title(title, fontsize = 20)
    fig.tight_layout()
    plt.show()
<<<<<<< Updated upstream
    plt.axis('off')
    if save:
        fig.savefig(os.path.join(filepath, title + \
                                 '.png'), dpi=300)
=======
    if save:
        fig.savefig(os.path.join(results_folder, title +'.png'))
>>>>>>> Stashed changes
    return

def plot_efficiency_percentile(efficiency, colors, results_folder = '', save = False):
    fig, ax = plt.subplots(2, sharex = True)
    for e in efficiency.keys():
        ax[0].plot(efficiency[e]['global_efficiency'][1], efficiency[e]['global_efficiency'][0], color = colors[e])
        ax[0].set_title('Global Efficiency')
        ax[1].plot(efficiency[e]['local_efficiency'][1], efficiency[e]['local_efficiency'][0], color = colors[e])
        ax[1].set_title('Local Efficiency')
        ax[1].set_xlabel('Percentile')
        start, end = ax[1].get_xlim()
        ax[1].xaxis.set_ticks(np.arange(start, end, 4.1))
    fig.legend(handles = create_custom_legend(colors))
    if save:
        fig.savefig(os.path.join(results_folder, 'Efficiency per percentile' + \
                                 '.png'), dpi=300)
    return None
def color_allens(allens_unique):
    color_list = [color for color in sns.color_palette('Set3', n_colors=len(allens_unique))]
    color_dict = {r: color_list[i] for i, r in enumerate(allens_unique)}
    return color_dict
def create_custom_legend(colors):
    """
    Creates a personalized legend taking into account a dictionary with the strings and the colors

    Parameters
    ----------
    colors : dict
        Dictionary in which the key is the string you want to use in the legend and the value is the color you want to assing.

<<<<<<< Updated upstream
=======
def plot_network_statistic(a):
    return

def plot_efficiency_percentile(efficiency, colors_exp, results_folder = '', save = False):
    fig, ax = plt.subplots(2, figsize = (20,16),sharex = True)
    for e in efficiency.keys():
        ax[0].plot(efficiency[e]['global_efficiency'][1], efficiency[e]['global_efficiency'][0], color = colors_exp[e],linewidth=3)
        ax[0].set_title('Global Efficiency')
        ax[1].plot(efficiency[e]['local_efficiency'][1], efficiency[e]['local_efficiency'][0], color = colors_exp[e],linewidth=3)
        ax[1].set_title('Local Efficiency')
        ax[1].set_xlabel('Percentile')
        #start, end = ax[1].get_xlim()
        #ax[1].xaxis.set_ticks(np.arange(start, end, 4.1))
    fig.legend(handles = create_custom_legend(colors_exp), fontsize = 'x-small')
    if save:
        fig.savefig(os.path.join(results_folder, 'Efficiency per percentile' + \
                                 '.png'), dpi=300)
    return None

def plot_modularity(modularity_dict, colors, results_folder = '', save = False):
    fig = plt.figure()
    for e,m in modularity_dict.items():
        plt.plot(m[0], m[1], color = colors[e])
        plt.title('Modularity')
        plt.xlabel('Percentile')
    fig.legend(handles = create_custom_legend(colors))
    if save:
        fig.savefig(os.path.join(results_folder, 'Modularity per percentile' + \
                                 '.png'), dpi=300)
    return None


def color_allens(allens_unique):
    color_list = [color for color in sns.color_palette('Set3', n_colors=len(allens_unique))]
    color_dict = {r: color_list[i] for i, r in enumerate(allens_unique)}
    return color_dict
def create_custom_legend(colors):
    """
    Creates a personalized legend taking into account a dictionary with the strings and the colors

    Parameters
    ----------
    colors : dict
        Dictionary in which the key is the string you want to use in the legend and the value is the color you want to assing.

>>>>>>> Stashed changes
    Returns
    -------
    handles
        Handles elements prepared to be plotted with matplotlib.legend.

    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    legend_elements = []
    for area, color in colors.items():
        legend_elements.append(Line2D([0], [0], color=color, lw=4, label = area))
    return legend_elements