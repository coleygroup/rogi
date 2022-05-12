#!/usr/bin/env python
import numpy as np
from copy import deepcopy
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import squareform, pdist
from numpy.typing import ArrayLike, NDArray


def MODI(Dx: NDArray, Y: ArrayLike) -> float:
    """
    Modelability index (MODI) for binary classification tasks, as described in:
    Golbraikh et al. "Data Set Modelability by QSAR", J. Chem. Inf. Model. 2014, 54, 1−4
    (https://pubs.acs.org/doi/10.1021/ci400572x).

    Parameters
    ----------
    Dx : ndarray
        Distance matrix in square form.
    Y : array_like
        Binary labels, 0 or 1.

    Returns
    -------
    float
        MODI score.
    """

    Y = np.array(Y)

    # make diagonal largest value, i.e. no self neighbors
    Dx[np.diag_indices(Dx.shape[0])] = np.max(Dx)

    # get index of neighbors
    neigh = NearestNeighbors(n_neighbors=1, metric='precomputed')
    neigh.fit(Dx)
    neighbors = neigh.kneighbors(X=Dx, n_neighbors=1, return_distance=False)
    neigh_indices = neighbors.flatten()  # idx -> idx of neighbor

    # get indices of elements of the two classes
    class_0_indices = np.where(Y < 0.5)[0]
    class_1_indices = np.where(Y > 0.5)[0]

    def get_fraction_of_neighbors_in_same_class(indices, neighbor_map):
        # for all compounds of a class:
        #   get nearest neighbor
        #   if same class -> +1
        # return n / |class|
        n = 0
        for i in indices:
            neigh_idx = neighbor_map[i]
            if neigh_idx in indices:
                n += 1
        return n / len(indices)

    F0 = get_fraction_of_neighbors_in_same_class(class_0_indices, neigh_indices)
    F1 = get_fraction_of_neighbors_in_same_class(class_1_indices, neigh_indices)

    modi_score = 0.5 * (F0 + F1)
    return modi_score


def RMODI(Dx: NDArray, Y: ArrayLike, delta: float = 0.625) -> float:
    """
    Regression Modelability index (RMODI) for regression, as described in:
    Ruiz et al. 2018. “Regression Modelability Index: A New Index for Prediction of the Modelability of Data Sets in
    the Development of QSAR Regression Models.” Journal of Chemical Information and Modeling 58 (10): 2069–84.
    (https://pubs.acs.org/doi/10.1021/acs.jcim.8b00313)

    Parameters
    ----------
    Dx : ndarray
        Distance matrix in square form.
    Y : array_like
        Continuous property values.
    delta : float
        Delta value used for Eq.9 in Ruiz et al. 2018. Default is 0.625.

    Returns
    -------
    float
        RMODI score.
    """

    # make diagonal largest value, i.e. no self neighbors
    Dx[np.diag_indices(Dx.shape[0])] = np.max(Dx) + 0.1

    # get property distance matrix
    Y = np.array(Y)
    y_std = Y.std()
    if Y.ndim < 2:
        Y = Y.reshape(-1, 1)
    Dy = squareform(pdist(Y, metric='euclidean'))

    # indices of mols in same/different classes
    idx_same = np.where(Dy <= delta * y_std)
    idx_diff = np.where(Dy > delta * y_std)

    # mask distance matrices
    Dx_same = deepcopy(Dx)
    Dx_same[idx_diff] = np.max(Dx) + 0.1
    Dx_diff = deepcopy(Dx)
    Dx_diff[idx_same] = np.max(Dx) + 0.1

    # get index of neighbors
    neigh = NearestNeighbors(n_neighbors=1, metric='precomputed')
    neigh.fit(Dx_same)
    neighbors = neigh.kneighbors(X=Dx_same, n_neighbors=1, return_distance=False)
    idx2neigh_same = neighbors.flatten()  # idx -> idx of neighbor in same class

    neigh = NearestNeighbors(n_neighbors=1, metric='precomputed')
    neigh.fit(Dx_diff)
    neighbors = neigh.kneighbors(X=Dx_diff, n_neighbors=1, return_distance=False)
    idx2neigh_diff = neighbors.flatten()  # idx -> idx of neighbor in different class

    # we have all the info
    # for each molecule, compute rivality index and add 1 if it is less than zero
    # i.e. if the closest molecule is in same class
    rmodi = 0.
    for i in range(Dx.shape[0]):
        idx_same = idx2neigh_same[i]
        idx_diff = idx2neigh_diff[i]
        dist_same = Dx[i, idx_same]
        dist_diff = Dx[i, idx_diff]
        rivality_index = (dist_same - dist_diff) / (dist_same + dist_diff)
        if rivality_index < 0:
            rmodi += 1.

    return rmodi / Dx.shape[0]