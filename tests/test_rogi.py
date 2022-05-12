#!/usr/bin/env python
from rogi import RoughnessIndex
import math
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np
from scipy.spatial.distance import squareform


def test_smiles_input():
    Y = [0, 1, 2, 3, 4]
    smiles = ['C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1',
              'C[C@H]1C[C@H]1C(=O)N1CCN(C(=O)NCC(C)(C)[NH+]2CCCCC2)CC1',
              'Cc1ccc([C@@H]2CSCCN2Cc2cc3cnn(C(C)C)c3nc2Cl)o1',
              'C=CCc1ccccc1OC(C)=O',
              'C[C@@H](c1nccs1)N(C)C(=O)c1cccc(NC(=O)C2CCCC2)c1']
    ri = RoughnessIndex(Y=Y, smiles=smiles)
    score = ri.compute_index()
    assert math.isclose(score, 0.2722166401201682)


def test_fingerprints():
    Y = [0, 1, 2, 3, 4]
    smiles = ['C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1',
              'C[C@H]1C[C@H]1C(=O)N1CCN(C(=O)NCC(C)(C)[NH+]2CCCCC2)CC1',
              'Cc1ccc([C@@H]2CSCCN2Cc2cc3cnn(C(C)C)c3nc2Cl)o1',
              'C=CCc1ccccc1OC(C)=O',
              'C[C@@H](c1nccs1)N(C)C(=O)c1cccc(NC(=O)C2CCCC2)c1']
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=2048) for m in mols]
    ri = RoughnessIndex(Y=Y, fps=fps)
    score = ri.compute_index()
    assert math.isclose(score, 0.2722166401201682)


def test_descriptors():
    Y = [0, 1, 2, 3, 4]
    smiles = ['C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1',
              'C[C@H]1C[C@H]1C(=O)N1CCN(C(=O)NCC(C)(C)[NH+]2CCCCC2)CC1',
              'Cc1ccc([C@@H]2CSCCN2Cc2cc3cnn(C(C)C)c3nc2Cl)o1',
              'C=CCc1ccccc1OC(C)=O',
              'C[C@@H](c1nccs1)N(C)C(=O)c1cccc(NC(=O)C2CCCC2)c1']

    descriptor_funcs = {desc[0]: desc[1] for desc in Descriptors.descList}
    X = []
    for smi in smiles:
        m = Chem.MolFromSmiles(smi)
        Xi = []
        for desc_name in ['MolLogP', 'qed', 'TPSA']:
            desc = descriptor_funcs[desc_name](m)
            Xi.append(desc)
        X.append(Xi)
    X = np.array(X)
    X = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    ri = RoughnessIndex(Y=Y, X=X, metric='euclidean')
    score = ri.compute_index()
    assert math.isclose(score, 0.3068390625203159)


def test_distance_matrix():
    Y = [0, 1, 2, 3, 4]
    d = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]) / 10.
    X = squareform(d)
    ri = RoughnessIndex(Y=Y, X=X, metric='precomputed')
    score = ri.compute_index()
    assert math.isclose(score, 0.1530913965143058)
