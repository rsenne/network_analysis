# import necessary libraries
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.special as sc
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.patches as mpatches
from tqdm import tqdm
import igraph as ig
from bct.algorithms import centrality
from scipy.spatial.distance import cdist

__all__ = ["NetworkAnalysis", "compare_spectrum", "plot_and_compare_degree_distribution"]


class NetworkAnalysis:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data = None
        self.G = None
        self.nodes = None
        self.node_attrs_df = None

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        self.data = self.data.apply(lambda x: x.fillna(x.mean()), axis=0)
        node_names = self.data.columns.to_list()
        self.nodes = {i: name for i, name in enumerate(node_names)}

    def corr_matrix(self, data, corr_type='Pearson', z_trans=True):
        if corr_type == 'Pearson':
            r_val, p_val = self.corr_pearson(data)
            p_adjusted = self.apply_correction(p_val)
            if z_trans:
                return np.arctanh(r_val), p_adjusted
            else:
                return r_val, p_adjusted
        elif corr_type == 'Spearman':
            r_val, p_val = stats.spearmanr(data, axis=0)
            p_adjusted = self.apply_correction(p_val)
            return r_val, p_adjusted

    def corr_pearson(self, data):
        dim = data.shape[1]
        r_val = np.zeros((dim, dim))
        p_val = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(i + 1, dim):
                r, p = stats.pearsonr(data[:, i], data[:, j])
                r_val[i, j] = r_val[j, i] = r
                p_val[i, j] = p_val[j, i] = p
        return r_val, p_val

    def apply_correction(self, p_val):
        mask = np.triu(np.ones(p_val.shape), 1).astype(
            bool)  # Mask to get upper triangular matrix excluding diagonal
        p_val_half = p_val[mask]
        _, p_adjusted_half, _, _ = multipletests(p_val_half, alpha=0.05, method='fdr_bh')
        p_adjusted = np.zeros_like(p_val)
        p_adjusted[mask] = p_adjusted_half
        p_adjusted = p_adjusted + p_adjusted.T - np.diag(np.diag(p_adjusted))  # Reconstruct full matrix
        return p_adjusted

    def percentile(self, array, p):
        num_obs = int(np.size(array, 0) ** 2 * p)
        crit_value = -np.sort(-array.flatten())[num_obs - 1]
        percent_arr = np.where(array < crit_value, 0, array)
        return percent_arr

    def significance_check(self, p_adjusted, corr, alpha, threshold=0.0, include_negs=True):
        p_adjusted = np.where((p_adjusted >= alpha), 0, p_adjusted)
        np.fill_diagonal(p_adjusted, 0)
        p_adjusted = np.where((p_adjusted != 0), 1, p_adjusted)
        if not include_negs:
            p_adjusted = np.where((corr < 0), 0, p_adjusted)
        threshold_matrix = np.multiply(corr, p_adjusted)
        threshold_matrix = np.where((abs(threshold_matrix) < threshold), 0, threshold_matrix)
        return threshold_matrix

    def networkx(self, corr_data, node_label, drop_islands=False):
        graph = nx.from_numpy_array(corr_data, create_using=nx.Graph)
        if node_label:
            graph = nx.relabel_nodes(graph, node_label)
        if drop_islands:
            remove_node = [node for node, degree in graph.degree() if degree < 1]
            graph.remove_nodes_from(remove_node)
        return graph

    def create_network(self, corr_type='Spearman', p_threshold=0.001, global_threshold=0.25, drop_islands=False):
        self.load_data()
        rvals, p = self.corr_matrix(self.data, corr_type=corr_type, z_trans=False)
        threshold_matrix = self.significance_check(p, rvals, alpha=p_threshold, include_negs=True)
        per_matrix = self.percentile(threshold_matrix, p=global_threshold)
        self.G = self.networkx(threshold_matrix, self.nodes, drop_islands=drop_islands)

    def compute_spectrum(self, G=None):
        G = G or self.G
        A = nx.adjacency_matrix(G)
        return np.linalg.eigvals(A.todense())

    def get_node_attributes(self, graph, use_distance=False, compress_to_df=False):
        if use_distance:
            G_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in graph.edges(data='weight')}
            nx.set_edge_attributes(graph, values=G_distance_dict, name='distance')
        deg = nx.degree_centrality(graph)
        between = nx.betweenness_centrality(graph)
        eig = nx.eigenvector_centrality(graph)
        close = nx.closeness_centrality(graph)
        clust = nx.clustering(graph)
        deg_sort = {area: val for area, val in sorted(deg.items(), key=lambda ele: ele[0])}
        between_sort = {area: val for area, val in sorted(between.items(), key=lambda ele: ele[0])}
        eig_sort = {area: val for area, val in sorted(eig.items(), key=lambda ele: ele[0])}
        close_sort = {area: val for area, val in sorted(close.items(), key=lambda ele: ele[0])}
        clust_sort = {area: val for area, val in sorted(clust.items(), key=lambda ele: ele[0])}
        if compress_to_df:
            node_info = {
                'Degree': list(deg_sort.values()),
                'Betweenness': list(between_sort.values()),
                'Eigenvector_Centrality': list(eig_sort.values()),
                'Closeness': list(close_sort.values()),
                'Clustering_Coefficient': list(clust_sort.values()),
            }
            ROI_index = list(graph.nodes)
            node_attrs_df = pd.DataFrame(node_info, index=ROI_index, columns=node_info.keys())
            return node_attrs_df
        else:
            return deg_sort, between_sort, eig_sort, close_sort, clust_sort

    def find_hubs(self, node_attrs_df):
        Results = node_attrs_df
        Results['Hub_Score'] = 0
        Results['Hub_Score'] = np.where((Results['Degree'] >= Results.Degree.quantile(0.80)), Results['Hub_Score'] + 1,
                                        Results['Hub_Score'])
        Results['Hub_Score'] = np.where(
            (Results['Eigenvector_Centrality'] >= Results.Eigenvector_Centrality.quantile(.80)),
            Results['Hub_Score'] + 1,
            Results['Hub_Score'])
        Results['Hub_Score'] = np.where((Results['Betweenness'] >= Results.Betweenness.quantile(0.80)),
                                        Results['Hub_Score'] + 1,
                                        Results['Hub_Score'])
        Results['Hub_Score'] = np.where(
            (Results['Clustering_Coefficient'] <= Results.Clustering_Coefficient.quantile(.20)),
            Results['Hub_Score'] + 1,
            Results['Hub_Score'])
        Results['Hub_Score'] = np.where((Results['Closeness'] >= Results.Closeness.quantile(.80)),
                                        Results['Hub_Score'] + 1, Results['Hub_Score'])

        NonHubs = Results[
            (Results['Hub_Score'] < 3)].index

        Hubs = Results.drop(NonHubs).sort_values('Hub_Score', ascending=False)
        return Results, Hubs

    def threshold_simulation(self, adj_mat, a, b, x, algo='markov'):
        percentiles = [i for i in np.linspace(a, b, x)]
        thresholded_arrays = [self.percentile(adj_mat, p) for p in percentiles]
        if algo == 'markov':
            modularity = []
            for thresh in thresholded_arrays:
                G, _ = self.networkx(thresh, node_label=None)
                _, mc_clusters = nx.algorithms.community.markov(G)
                modularity.append(nx.algorithms.community.modularity(G, mc_clusters))
        else:
            modularity = []
        return percentiles, modularity

    def combine_node_attrs(self, node_attrs_df, WMDz_PC_df, Allens):
        final_df = pd.merge(node_attrs_df, WMDz_PC_df, left_index=True, right_index=True)
        final_df['Allen_ROI'] = Allens
        final_df = final_df[
            ["Allen_ROI", "Degree", "Betweenness", "Eigenvector_Centrality", "Closeness", "Clustering_Coefficient",
             "WMDz",
             "PC", "Hub_Score"]]
        return final_df

    def leiden_algorithm(self, resolution_range=(0.5, 1.75), resolution_steps=1000, n_iterations=1000):
        if self.G is None:
            raise ValueError("Please create the network graph using create_network method first.")

        # Convert networkx graph to igraph
        ig_graph = ig.Graph.from_networkx(self.G)

        mod = []
        resolutions = np.linspace(*resolution_range, resolution_steps)
        for i in tqdm(resolutions):
            partition = ig_graph.community_leiden(objective_function='CPM', weights='weight', resolution=i,
                                                  n_iterations=n_iterations)
            mod.append(partition.graph.modularity(partition.membership))
        max_mod_idx = np.argmax(mod)
        best_resolution = resolutions[max_mod_idx]

        best_partition = ig_graph.community_leiden(objective_function='CPM', weights='weight',
                                                   resolution=best_resolution, n_iterations=n_iterations)

        grouped_indices = {value: [i for i, x in enumerate(best_partition.membership) if x == value] for value in
                           set(best_partition.membership)}
        community_assignment = [tuple(indices) for indices in grouped_indices.values()]

        return best_partition, community_assignment, best_resolution


# Will generate a euclidean distance matrix from the raw data
def eucl_matrix(data):
    data = data.T
    eucl_matrix = cdist(data, data, metric='euclidean')
    return eucl_matrix


def shortest(G, threshold_matrix):
    # Function to calculate shortest path of each node
    short = nx.floyd_warshall_numpy(G, weight='weight')
    shortavg = np.mean(short, axis=0)
    keys = threshold_matrix.index.values.tolist()
    vals = shortavg.tolist()
    zip_iterator = zip(keys, vals)
    short_dictionary = dict(zip_iterator)

    return short_dictionary


def get_ordered_list(G, stat='Degree'):
    if stat == 'Degree':
        ordered = {k: v for k, v in sorted(dict(G.degree()).items(), key=lambda item: item[1])}
    elif stat == 'WeightedDegree':
        ordered = {k: v for k, v in sorted(dict(G.degree(weight='weight')).items(), key=lambda item: item[1])}
    elif stat == 'EigCentrality':
        ordered = {k: v for k, v in sorted(dict(nx.eigenvector_centrality(G)).items(), key=lambda item: item[1])}
    elif stat == 'Clustering':
        ordered = {k: v for k, v in sorted(dict(nx.clustering(G)).items(), key=lambda item: item[1])}
    elif stat == 'Betweenness':
        ordered = {k: v for k, v in sorted(dict(nx.betweenness_centrality(G)).items(), key=lambda item: item[1])}
    else:
        ordered = []
    return ordered


def cluster_attributes(graph, nodes, communities, make_df=False):
    adj_matrix = nx.to_numpy_array(graph)  # Will create an adjacency matrix from the graph as a np.ndarray
    node_ROIs = nodes.values()
    WMDz = centrality.module_degree_zscore(adj_matrix, communities, flag=0)  # calculate the WMDz
    PC = centrality.participation_coef(adj_matrix, communities, 'undirected')  # calculate the participation coefficient
    if make_df:
        d = {'WMDz': WMDz, "PC": PC}
        df = pd.DataFrame(d, columns=["WMDz", "PC"], index=node_ROIs)
        return df
    else:
        return WMDz, PC


def node_attrs_to_csv(final_df, folder, var_name):
    final_df.to_csv(folder + '/' + var_name + '.csv')
    return


def plot_and_compare_degree_distribution(G1, G2):
    degree_sequence1 = sorted([d for n, d in G1.degree()], reverse=True)
    degree_sequence2 = sorted([d for n, d in G2.degree()], reverse=True)

    plt.figure(figsize=(10, 8))
    plt.plot(degree_sequence1, 'b-', marker='o', label="G1")
    plt.plot(degree_sequence2, 'r-', marker='o', label="G2")
    plt.title("Degree Rank Plot")
    plt.ylabel("Degree")
    plt.xlabel("Rank")
    plt.legend()
    plt.show()

    # KS test
    ks_stat, p_value = stats.ks_2samp(degree_sequence1, degree_sequence2)
    print(f"KS statistic: {ks_stat}")
    print(f"p-value: {p_value}")


def compare_spectrum(G1, G2):
    spectrum1 = G1.compute_spectrum()
    spectrum2 = G2.compute_spectrum()
    return np.linalg.norm(spectrum1 - spectrum2)


def corr_permutation_ttest(a, b, niters=1000, plot=False):
    upper_idxs = np.triu_indices_from(a, 1)
    upper_a = a[upper_idxs]
    upper_b = b[upper_idxs]
    mean_a = np.mean(upper_a)
    mean_b = np.mean(upper_b)
    test_statistic = np.abs(mean_a - mean_b)

    combined = np.concatenate([upper_a, upper_b])
    num_a = len(upper_a)
    num_b = len(upper_b)
    perm_diffs = np.zeros(niters)

    for i in range(niters):
        permuted = np.random.permutation(combined)
        perm_a = permuted[:num_a]
        perm_b = permuted[num_a:num_a + num_b]
        perm_diffs[i] = np.abs(np.mean(perm_a) - np.mean(perm_b))

    # p-value is proportion of permuted differences greater than or equal to observed difference
    p_value = np.mean(perm_diffs >= test_statistic)

    # Plotting
    if plot:
        plt.hist(perm_diffs, bins=30, alpha=0.7, label='Null distribution')
        plt.axvline(test_statistic, color='red', linestyle='dashed', linewidth=2, label='Observed statistic')
        plt.title("Null distribution vs observed statistic")
        plt.xlabel("Difference in means")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()

    return test_statistic, p_value
