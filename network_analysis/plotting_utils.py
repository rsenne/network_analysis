import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from matplotlib import pyplot as plt


def publication_heatmap(adj_mat, labels, save=True):
    fig, ax = plt.subplots()
    sns.heatmap(adj_mat, annot=labels, cmap='vlag', linewidths=0.3, ax=ax)
    plt.show()
    if save:
        plt.savefig()
    return


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
    allen_list = list(set(allen_df['Allen Area']))
    allen_list_alphabetical = sorted(allen_list)
    allen_colors = [color for color in sns.color_palette('Set3', len(allen_list_alphabetical))]
    allen_color_dict = {group: color for group, color in zip(allen_list_alphabetical, allen_colors)}
    allen_dict = {abbrev: name for abbrev, name in zip(allen_df['Abbreviation'], allen_df['Allen Area'])}
    color_list = []
    for area in allen_dict:
        color_list.append(allen_color_dict[allen_dict[area]])
    return color_list


def sunflower_theta(n):
    golden_ratio = ((1 + 5 ** 0.5) / 2) ** 2
    return 2 * np.pi / golden_ratio * n


def sunflower_r(n, c=1.8):
    return c * (n ** 0.5)


def get_point_cloud(k, c=1.8):
    n = [i for i in range(1, k + 1)]
    r = [sunflower_r(i, c) for i in n]
    theta = [sunflower_theta(j) for j in n]
    point_cloud_x = (r * np.cos(theta))
    point_cloud_y = (r * np.sin(theta))
    point_cloud = [[x, y] for x, y in zip(point_cloud_x, point_cloud_y)]
    return point_cloud


def get_position_data(cluster_list, node_names, shape='circular'):
    number_of_clusters = len(cluster_list)
    nodes_list = [x for x in range(0, number_of_clusters)]
    pos_graph = nx.Graph()
    pos_graph.add_nodes_from(nodes_list)
    if shape == 'circular':
        pos = nx.circular_layout(pos_graph, scale=39, dim=2)
    else:
        pos = {i: get_point_cloud(number_of_clusters, 12)[i] for i in range(number_of_clusters)}
    num_of_nodes = [len(node) for node in cluster_list]
    point_clouds = [get_point_cloud(num_of_nodes[lens], 1.8) for lens in range(len(num_of_nodes))]
    for i in range(len(point_clouds)):
        for j in range(len(point_clouds[i])):
            point_clouds[i][j][0] += pos[i][0]
            point_clouds[i][j][1] += pos[i][1]
    # avg_x_y = [[np.mean([point_cloud[i][0] for i in range(len(point_cloud))]), np.mean([point_cloud[i][1] for i in range(len(point_cloud))])]for point_cloud in point_clouds]
    # radii = [np.sqrt((i[0][0]-i[-1][0])**2 + (i[-1][1]-i[0][1])**2) for i in point_clouds]
    point_cloud_map = {cluster: pos_list for cluster, pos_list in enumerate(point_clouds)}
    pos_dict = {}
    for i in range(len(cluster_list)):
        for j in range(len(cluster_list[i])):
            pos_dict.update({node_names[cluster_list[i][j]]: np.array(point_cloud_map[i][j])})
    nx.rescale_layout_dict(pos_dict, 10)
    return pos_dict


def graph_network(G, color_list, pos_dict):
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 750 for degree in dict(deg).values()]
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos_dict, width=1, edge_color=edge_colors, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(G, pos=pos_dict, node_size=node_sizes, node_color=color_list, linewidths=1,
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos=pos_dict)
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.axis('off')
    plt.show()
    return fig


def plot_degree_distribution(G):
    fig, ax = plt.subplots()
    degree_values = list(dict(G.degree()).values())
    kde = stats.gaussian_kde(degree_values)
    kde_lin = np.linspace(np.min(degree_values), np.max(degree_values), 1000)
    ax.hist(degree_values, density=True)
    ax.plot(kde_lin, kde(kde_lin))
    ax.plot()
    return


def plot_r_distributions(adj1, adj2):
    fig, ax = plt.subplots()
    bins = np.arange(-1, 7, 0.1)
    ax.hist(adj1.flatten(), density=True, bins=bins, align='mid', color='cornflowerblue', alpha=0.7)
    ax.hist(adj2.flatten(), density=True, bins=bins, align='mid', color='tomato', alpha=0.7)
    kde = stats.gaussian_kde(adj1.flatten())
    kde2 = stats.gaussian_kde(adj2.flatten())
    kde_lin = np.linspace(-1, 6, 1000)
    ax.plot(kde_lin, kde(kde_lin), color='cornflowerblue')
    ax.plot(kde_lin, kde2(kde_lin), color='tomato')
    return


'''def plot_network_statistic(G, ):
    return'''
