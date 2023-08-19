import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import math
from scipy import stats
import random
from matplotlib import pyplot as plt


def publication_heatmap(adj_mat, labels, save=True):
    fig, ax = plt.subplots()
    sns.heatmap(adj_mat, cmap='vlag', linewidths=0.3, ax=ax)
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


def sunflower_r(n, c=2.0):
    return c * (n ** 0.5)


def get_point_cloud(k, c=2.0):
    n = np.arange(1, k + 1)
    r = sunflower_r(n, c)
    theta = sunflower_theta(n)
    point_cloud_x = r * np.cos(theta)
    point_cloud_y = r * np.sin(theta)
    return np.column_stack((point_cloud_x, point_cloud_y))


def generate_circle_centers(num_circles, min_dist, radii, max_attempts=1000, max_space=500):
    """Generate random (x, y) coordinates for circle centers with given radii"""
    space = 175

    while space <= max_space:
        centers = [(random.uniform(0, space), random.uniform(0, space))]

        attempts = 0
        while len(centers) < num_circles and attempts < max_attempts:
            # Generate a new random center
            new_center = (random.uniform(0, space), random.uniform(0, space))

            # Check if the new circle overlaps with any existing circle
            too_close = any(np.hypot(*(np.array(center) - np.array(new_center))) < min_dist + radius
                            for center, radius in zip(centers, radii))

            # Add the new circle if it doesn't overlap
            if not too_close:
                centers.append(new_center)
            attempts += 1

        if len(centers) == num_circles:
            return centers  # Successfully placed all circles, return the centers

        print(
            f"Warning: could only place {len(centers)} out of {num_circles} circles within {max_attempts} attempts. Increasing space.")
        space += 10  # Increase space by 10 units for the next iteration

    # If we reach here, it means we couldn't place all circles even after increasing the space
    print(
        f"Error: could not place all circles even with space={space - 10}. Placed only {len(centers)} out of {num_circles}.")
    return centers


def get_position_data(cluster_list, node_names, shape='circular'):
    number_of_clusters = len(cluster_list)
    num_of_nodes = [len(node) for node in cluster_list]
    pos = generate_circle_centers(number_of_clusters, 40, np.array(num_of_nodes) * (1 / 10))
    point_clouds = [get_point_cloud(num, 7.5) + np.array(center) if num > 10 else get_point_cloud(num, 7) + np.array(center) for num, center in zip(num_of_nodes, pos)]

    pos_dict = {node_names[node]: point for idx, cluster in enumerate(cluster_list) for node, point in
                zip(cluster, point_clouds[idx])}

    return pos_dict

def graph_network(G, color_list, cluster_list, node_names, figsize=(25, 25)):
    # Filter out communities with nodes of degree 0
    cluster_list = [cluster for cluster in cluster_list if all(G.degree(node_names[node]) != 0 for node in cluster)]
    print(cluster_list)
    
    # Filter node_names based on the filtered cluster_list
    filtered_node_names = {node: node_names[node] for cluster in cluster_list for node in cluster if node in node_names}
    
    # Generate position data for the remaining clusters
    pos_dict = get_position_data(cluster_list, filtered_node_names)

    fig, ax = plt.subplots(figsize=figsize)
    negativeCorr, positiveCorr = 'lightcoral', 'gainsboro'
    
    # Filter edges based on nodes in pos_dict
    filtered_edges = [(i, j) for i, j in G.edges if i in pos_dict and j in pos_dict]
    edge_colors = [negativeCorr if G[i][j]['weight'] < 0 else positiveCorr for i, j in filtered_edges]

    # Filter nodes and compute node sizes
    filtered_nodes = [node for node in G.nodes if node in pos_dict]

    # Filter colors
    reverse_node_names = {v: k for k, v in node_names.items()}
    filtered_node_indices = [reverse_node_names[node] for node in filtered_nodes]
    filtered_colors = [color_list[idx] for idx in filtered_node_indices]
    # node sizes
    node_sizes = [G.degree(node) / np.mean([G.degree(n) for n in filtered_nodes]) * 500 for node in filtered_nodes]
    H = G.subgraph(filtered_nodes)
    nx.draw(H, pos=pos_dict,
            nodelist=filtered_nodes,
            edgelist=filtered_edges,
            node_color=filtered_colors,
            node_size=node_sizes,
            linewidths=1,
            edge_color=edge_colors,
            edgecolors='black',
            with_labels=True,
            ax=ax)
    plt.axis('off')
    return fig


def DG_subgraph(cluster_ids, nodes, G, pos_dict, color_list):
    DG_tuple = [clust for clust in cluster_ids if 24 in clust]
    DG_cluster = [node for tup in DG_tuple for node in tup]

    vals = []
    for node in DG_cluster:
        vals.append(nodes[node])

    new_color_list = []
    for node in DG_cluster:
        new_color_list.append(color_list[node])

    sub_graph = G.subgraph(vals)
    edge_list = list(sorted(sub_graph.edges()))

    DG_graph = nx.Graph(pos=pos_dict)
    DG_graph.add_nodes_from(sorted(G.subgraph(vals)))
    DG_graph.add_edges_from((edge_list))

    deg = DG_graph.degree()
    node_sizes = [degree / np.mean(list(dict(deg).values())) * 1000 for degree in dict(deg).values()]
    fig, ax = plt.subplots(figsize=(20, 15))
    nx.draw_networkx_edges(DG_graph, pos=pos_dict, width=1, edge_color='gainsboro', connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_nodes(DG_graph, pos=pos_dict, node_size=node_sizes, node_color=new_color_list, linewidths=1,
                           edgecolors='black')
    nx.draw_networkx_labels(DG_graph, pos=pos_dict)
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.axis('off')
    plt.show()
    return DG_graph


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


def plot_network_statistic(node_attr_df, statistic: str = 'Degree', cutoff: float = 0.2):
    stat = node_attr_df.loc[node_attr_df[statistic] < node_attr_df[statistic].quantile(cutoff)][statistic].sort_values(
        ascending=True)
    rois = list(stat.index)
    fig, ax = plt.subplots(figsize=(6.4 * 3, 4.8))
    ax.bar(rois, stat, color='b', linewidth=0.3, edgecolor='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis="x", rotation=90)
    return fig, ax
