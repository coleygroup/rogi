#!/usr/bin/env python
from rogi import SARI
import math
import numpy as np
from scipy.spatial.distance import squareform


def test_smiles_input():
    Y = [6.1, 6.2, 6.3, 6.4, 6.5]
    smiles = ['c1ccccc1',
              'c1cc(C)ccc1',
              'Cc1ccc([C@@H]2CSCCN2Cc2cc3cnn(C(C)C)c3nc2Cl)o1',
              'C=CCc1ccccc1OC(C)=O',
              'C[C@@H](c1nccs1)N(C)C(=O)c1cccc(NC(=O)C2CCCC2)c1']
    sari = SARI(pKi=Y, smiles=smiles)
    score = sari.compute_sari()
    assert math.isclose(score, 0.1846109751180401)

    score = sari.compute_sari(mean_raw_cont=0.5, std_raw_cont=0.1, mean_raw_disc=-0.5, std_raw_disc=2.1)
    assert math.isclose(score, 0.19569841146521502)


def test_precomputed_distance_matrix():
    Y = [6.1, 6.2, 6.3, 6.4, 6.5]
    d = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]) / 10.
    sim_matrix = squareform(d)
    sari = SARI(pKi=Y, sim_matrix=sim_matrix)
    score = sari.compute_sari()
    assert math.isclose(score, 0.14831257148308755)
