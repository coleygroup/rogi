#!/usr/bin/env python
import numpy as np
from scipy import stats
from typing import Tuple, Optional
from numpy.typing import ArrayLike, NDArray

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, MACCSkeys, DataStructs
except ImportError:
    raise ImportError("cannot import rdkit, please install it via conda:\nconda install rdkit>=2021")


class SARI:
    def __init__(self, pKi: ArrayLike, smiles: Optional[ArrayLike] = None,
                 fingerprints: str = 'maccs', sim_matrix: Optional[NDArray] = None):
        """
        SAR Index as described in Peltason and Bajorath 2007, “SAR Index: Quantifying the Nature of Structure−Activity
        Relationships.” Journal of Medicinal Chemistry, 50 (23): 5571–78
        (https://pubs.acs.org/doi/pdf/10.1021/jm0705713).

        Parameters
        ----------
        pKi : array_like
            Array of log-affinity values, like pKi or pIC50.
        smiles : smiles
            Array of SMILES in the order that matches that of the affinity values provided.
        fingerprints : str
            Type of fingerprints to use. Available options are "morgan" or "maccs". Default is "maccs". Morgan
            fingerprints use length 2048 and radius 2.
        sim_matrix : array_like, optional
            Precomputed similarity matrix. If provided, the arguments `smiles` and `fingerprints` are ignored.
        """

        # co-domain distances
        self.y = pKi
        self.n_mols = len(pKi)

        # ------------------
        # domain distances
        # ------------------

        # if sim matrix is provided, we're done
        if sim_matrix is not None:
            self.sim_matrix = np.array(sim_matrix)
            assert self.sim_matrix.shape[0] == self.sim_matrix.shape[1]
            assert self.sim_matrix.shape[0] == self.n_mols
        # else, computer fingerptinrs and distances
        else:
            mols = [Chem.MolFromSmiles(smi) for smi in smiles]
            if isinstance(fingerprints, str):
                if fingerprints == 'maccs':
                    fps = [MACCSkeys.GenMACCSKeys(m) for m in mols]
                elif fingerprints == 'morgan':
                    fps = [AllChem.GetMorganFingerprintAsBitVect(m, nBits=2048, radius=2) for m in mols]
                else:
                    raise ValueError('only "maccs" or "morgan" fingerprints are allowed')
            else:
                fps = fingerprints

            assert len(fps) == self.n_mols

            # compute similarity matrix
            self.sim_matrix = np.zeros(shape=(self.n_mols, self.n_mols))
            for i in range(self.n_mols):
                # i+1 becauase we know the diagonal is zero
                self.sim_matrix[i, i + 1:] = np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
                self.sim_matrix[i + 1:, i] = self.sim_matrix[i, i + 1:]

    def compute_raw_scores(self) -> Tuple[float, float]:
        """Computes the raw SARI scores.

        Returns
        -------
        tuple
            The continuous and discontinuous raw scores.
        """
        # raw continuity score
        num = 0.
        den = 0.
        for i in range(self.n_mols):
            for j in range(self.n_mols):
                if i > j:
                    w_ij = (self.y[i] * self.y[j]) / (1 + np.abs(self.y[i] - self.y[j]))
                    den += w_ij
                    num += w_ij * self.sim_matrix[i, j]

        raw_cont = 1. - num / den

        # raw discontinuity score
        num = 0.
        den = 0.
        for i in range(self.n_mols):
            for j in range(self.n_mols):
                if i > j:
                    if self.sim_matrix[i, j] > 0.6:
                        num += np.abs(self.y[i] - self.y[j]) * self.sim_matrix[i, j]
                        den += 1.
        if den < 1e-6:
            raw_disc = 10 ** 8
        else:
            raw_disc = num / den

        return raw_cont, raw_disc

    def compute_sari(self, mean_raw_cont: float = 0., std_raw_cont: float = 1.,
                     mean_raw_disc: float = 0., std_raw_disc: float = 1.) -> float:
        """Computes the SARI score. The mean and standard deviation of raw scores across a reference set of
        datasets are used to standardize the SARI score, as described in https://pubs.acs.org/doi/pdf/10.1021/jm0705713.

        Parameters
        ----------
        mean_raw_cont : float
            Average of continuous raw scores. Default is zero.
        std_raw_cont : float
            Standard deviation of continuous raw scores. Default is one.
        mean_raw_disc : float
            Average of discontinuous raw scores. Default is zero.
        std_raw_disc : float
            Standard deviation of discontinuous raw scores. Default is one.

        Returns
        -------
        float
            SARI value.
        """
        raw_cont, raw_disc = self.compute_raw_scores()

        # z-scores
        zscore_cont = (raw_cont - mean_raw_cont) / std_raw_cont
        zscore_disc = (raw_disc - mean_raw_disc) / std_raw_disc

        # normalise by taking the CDF
        score_cont = stats.norm(loc=0.0, scale=1.0).cdf(zscore_cont)
        score_disc = stats.norm(loc=0.0, scale=1.0).cdf(zscore_disc)

        # SARI
        sari = (score_cont * (1. - score_disc)) * 0.5
        return sari