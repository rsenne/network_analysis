import networkx as nx
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt


def grab_color_attributes(cluster_list, node_dict):
    color_list = [color for color in sns.color_palette('colorblind', len(cluster_list))]
    color_dict = {}
    for i in range(0, len(cluster_list)):
        for j in cluster_list[i]:
            color_dict[node_dict[j]] = color_list[i]
    color_dict_sorted = {area: color for area, color in sorted(color_dict.items(), key=lambda ele: ele[0])}
    return color_dict_sorted


def sunflower_theta(n):
    golden_ratio = ((1 + 5 ** 0.5) / 2) ** 2
    return 2 * np.pi / golden_ratio * n


def sunflower_r(n, c=1):
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
    pos = nx.circular_layout(pos_graph, scale=30, dim=2)
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


def graph_network(G, color_list, pos_dict):
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j, _ in G.edges(data=True)]
    deg = G.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    # pos = nx.spring_layout(G, k=0.5, seed=3847897236)
    # print(pos)
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(G, pos=pos_dict, width=1, edge_color=edge_colors, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(G, pos=pos_dict, node_size=node_sizes, node_color=color_list, linewidths=1,
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos=pos_dict)
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.show()
    plt.axis('off')
    return
