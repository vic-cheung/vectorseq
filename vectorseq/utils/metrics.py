import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, euclidean
from itertools import combinations
from tqdm.auto import tqdm
from typing import Union, List


def sum_of_squares_error(numbers: Union[List[float], np.array]) -> float:
    """
    Compute the sum of squares error from a list of numbers.
    For a list of numbers [x_1, x_2, ... x_n]:
    sum of squares error = sum((x_i - mean) ** 2)
    """
    numbers = np.array(numbers)
    mean = np.mean(numbers)
    return np.sum((numbers - mean) ** 2)


def variance(numbers: Union[List[float], np.array]) -> float:
    """
    Compute variance from a list of numbers.
    For a list of numbers [x_1, x_2, ... x_n]:
    variance = sum((x_i - mean) ** 2) / (n - 1)
    """
    n = len(numbers)
    return sum_of_squares_error(numbers) / (n - 1)


def within_cluster_sum_of_squares(
    X: np.array, labels: Union[np.array, pd.Series]
) -> float:
    """
    Compute total within-cluster sum of squares error.

    For each cluster, the squared distance between each data point and the mean
    is computed.  We sum these distances to obtain within-cluster sum of square error.
    Then we add the within-cluster sum of square error across all clusters.

    Args:
        X: array-like of shape (n_samples, n_features).
        labels: array-like of shape (n_samples,).  Vector of cluster labels stored in adata.obs.

    Returns:
        Total sum of all clusters' within-cluster sum of squares error.  Sum of squares error
        is calculated for each of the clusters and then sum across all results is taken.
    """
    cluster_ids = sorted(list(set(labels)))
    cluster_ss_errs = []
    for cluster_id in cluster_ids:
        cluster_mask = labels == cluster_id
        cluster_X = X[cluster_mask, :]
        cluster_ss_errs += [sum_of_squares_error(cluster_X)]
    return np.sum(np.array(cluster_ss_errs))


def within_cluster_variance(X: np.array, labels: Union[np.array, pd.Series]) -> float:
    """
    Compute total within-cluster variance.

    Args:
        X: array-like of shape (n_samples, n_features).
        labels: array-like of shape (n_samples,).  Vector of cluster labels stored in adata.obs.

    Returns:
        Total sum of all clusters' within-cluster variance.  Sum of variance
        is calculated for each of the clusters and then sum across all results is taken.
    """
    cluster_ids = sorted(list(set(labels)))
    cluster_ss_errs = []
    for cluster_id in cluster_ids:
        cluster_mask = labels == cluster_id
        cluster_X = X[cluster_mask, :]
        cluster_ss_errs += [variance(cluster_X)]
    return np.sum(np.array(cluster_ss_errs))


def baker_hubert_gamma_index(X: np.array, labels: Union[np.array, pd.Series]) -> float:
    """
    Compute Baker-Hubert Gamma Index.


    Args:
        X: array-like of shape (n_samples, n_features).
        labels: array-like of shape (n_samples,).  Vector of cluster labels stored in adata.obs.

    Returns:
        Baker-Hubert Gamma Index

    References:
    1. Baker F and Hubert L. “Measuring the Power of Hierarchical Cluster Analysis”.
    Journal of the American Statistical Association 70 (1975): 31-38.
    2. Desgraupes, B. "Clustering Indices".  (2017)
    3. Vendramin L, Campello R, Hruschka E.  "Relative Clustering Validity Criteria:
    A Comparative Review". Statistical Analysis and Data Mining. (2010)
    4. Tomasini, C., Borges, E., Machado, K. and Emmendorfer, L. "A Study on the Relationship
    between Internal and External Validity Indices Applied to Partitioning and Density-based
    Clustering Algorithms". In Proceedings of the 19th International Conference on
    Enterprise Information Systems (2017)
    """
    # Generate Numeric Index for each of the Cluster Labels
    cluster_ids = np.unique(labels)
    num_clusters = len(cluster_ids)
    label2index = {k: v for k, v in zip(cluster_ids, range(num_clusters))}
    label_indices = np.array([label2index[label] for label in labels])

    # Proximity matrix, P(i,j) = pair-wise distance between observations i & j
    # Condensed representation is the concatenated rows of lower-triangular matrix.
    P_condensed = pdist(X, metric="euclidean")
    num_pairs = len(P_condensed)

    # Between vs. Within Cluster Membership
    # Y(i,j) = 0 if i,j are in same cluster, >0 if not in same cluster.
    temp = np.zeros((len(labels), 2))
    temp[:, 0] = label_indices
    Y_condensed = pdist(temp, metric="euclidean")

    # Compute s+ and s-
    # s+ = num times pair-wise dist within cluster < pair-wise dist outside cluster
    # s- = num times pair-wise dist within cluster > pair-wise dist outside cluster
    s_plus = 0
    s_minus = 0
    # Note: if pair-wise dist within cluster == pair-wise dist outside cluster, don't count it
    # Iterate through all unique pairs of pair-wise distance i,j in condensed Proximity matrix P

    # for i in tqdm(range(num_pairs - 1), desc="Pair i"):
    #     for j in tqdm(range(i + 1, num_pairs), desc="Pair j"):
    for i, j in tqdm(
        combinations(range(num_pairs), 2),
        total=num_pairs ** 2 - num_pairs,
        desc="Calculate s+,",
    ):
        # Distance i is within-cluster; Distance j is between-cluster
        if Y_condensed[i] == 0 and Y_condensed[j] > 0:
            within_clust_dist = P_condensed[i]
            between_clust_dist = P_condensed[j]
            if within_clust_dist < between_clust_dist:
                s_plus += 1
            if within_clust_dist > between_clust_dist:
                s_minus += 1
        # Distance i is between-cluster; Distance j is within-cluster
        if Y_condensed[i] > 0 and Y_condensed[j] == 0:
            within_clust_dist = P_condensed[j]
            between_clust_dist = P_condensed[i]
            if within_clust_dist < between_clust_dist:
                s_plus += 1
            if within_clust_dist > between_clust_dist:
                s_minus += 1

    # Compute Baker-Hubert's Gamma Index
    bh_gamma = (s_plus - s_minus) / (s_plus + s_minus)
    return bh_gamma


def xie_beni_index(X: np.array, labels: Union[np.array, pd.Series]) -> float:
    """
    Compute Xie-Beni Index.  A measure of separation vs. compactness.
    Smaller Xie-Beni Index indicates better cluster compactness.
    This implementation is for crisp (not fuzzy) clustering.

    XB = sum(|x_c - mean_c|^2) / (n * min_dist(c)^2)
    where
    x_c is a data point in a cluster c
    mean_c is the cluster center of cluster c
    n is the total number of data points
    min_dist(c) is the minimum distance between all cluster centers

    Args:
        X: array-like of shape (n_samples, n_features).
        labels: array-like of shape (n_samples,).  Vector of cluster labels stored in adata.obs.

    Returns:
        Xie-Beni Index.

    References:
    1. Xie X and Beni G. "A validity measure for fuzzy clustering".
    IEEE Transactions on Pattern Analysis and Machine Intelligence,
    (1991) 13:8, pp. 841-847
    2. Desgraupes B. "Clustering Indices".  (2017)
    3. Halkidi M, Batistakis Y and Vazirgiannis M. “On Clustering Validation
    Techniques.” Journal of Intelligent Information Systems 17 (2004): 107-145.
    """
    cluster_ids = list(set(labels))
    num_obs = len(labels)
    sum_mean_squared_norm = 0
    cluster_centers = []
    for cluster in cluster_ids:
        # Get data points in cluster
        indices = np.where(labels == cluster)[0]
        cluster_X = X[indices, :]
        # Cluster center
        cluster_center = np.mean(cluster_X, axis=0)
        cluster_centers += [cluster_center]
        # Mean squared euclidean distances of points in cluster w.r.t cluster center
        cluster_mean_squared_norms = [
            euclidean(x, cluster_center) ** 2 for x in cluster_X
        ]
        sum_mean_squared_norm += np.sum(cluster_mean_squared_norms)
    # Min distance between cluster centers
    centers = np.stack(cluster_centers)
    min_dist = min(pdist(centers))
    xb_index = sum_mean_squared_norm / (num_obs * (min_dist ** 2))
    return xb_index
