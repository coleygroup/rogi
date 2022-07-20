# Measures of roughness for molecular property landscapes

This package implements the roughness index (ROGI) presented in 
["Roughness of Molecular Property Landscapes and Its Impact on Modellability"](#), as well
as the [SARI](https://pubs.acs.org/doi/10.1021/jm0705713), [MODI](https://pubs.acs.org/doi/10.1021/ci400572x), 
and [RMODI](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00313) indices.

## Installation
``rogi`` can be installed with ``pip``:

```
pip install rogi
```

Note that ``rdkit`` is a dependency but needs to be installed separately with `conda`.

### Requirements
* `numpy`
* `scipy>=1.4`
* `fastcluster` 
* `pandas`
* `scikit-learn>=1`
* `rdkit >= 2021` to be installed with `conda`

## Usage
Note that ``ROGI`` and ``SARI`` are classes, while ``MODI`` and ``RMODI`` are functions. 

### ROGI
If SMILES are used as input, Morgan fingerprints (length 2048, radius 2) are computed and 
a distance matrix calculated with the Tanimoto metric:

```
from rogi import RoughnessIndex

ri = RoughnessIndex(Y=Y, smiles=smiles)
ri.compute_index()
>>> 0.42
```

With precomputed fingerprints:
```
ri = RoughnessIndex(Y=Y, fps=fingerprints)
ri.compute_index()
```

With descriptors you can pass a 2D array or a ``pandas.DataFrame`` where each row is a different
molecule, and each column a different descriptor:
```
ri = RoughnessIndex(Y=Y, X=descriptors, metric='euclidean')
ri.compute_index()
```

You can also precompute a distance matrix using any chosen representation and metric:
```
ri = RoughnessIndex(Y=Y, X=descriptors, metric='precomputed')
ri.compute_index()
```

### SARI
You can provide SMILES as input, and compute the SARI score without considering a reference
set of datasets as follows:
```
from rogi import SARI

sari = SARI(pKi=pKi, smiles=smiles, fingerprints='maccs')
sari.compute_sari()
>>> 0.42
```

To standardize the raw continuous and discontinuous scores based on a reference set of datasets,
you can compute the raw scores first and then provide SARI with their average and standard deviation:

```
raw_conts = []
raw_discs = []

for smiles, pKi in zip(datasets, affinities):
    sari = SARI(pKi=pKi, smiles=smiles, fingerprints='maccs')
    raw_cont, raw_disc = sari.compute_raw_scores()
    raw_conts.append(raw_cont)
    raw_discs.append(raw_disc)

mean_raw_cont = np.mean(raw_conts)
std_raw_cont = np.std(raw_conts)
mean_raw_disc = np.mean(raw_discs)
std_raw_disc = np.std(raw_discs)
                         
sari = SARI(pKi=my_pKi, smiles=my_smiles, fingerprints='maccs')
sari.compute_sari(mean_raw_cont=mean_raw_cont, std_raw_cont=std_raw_cont,
                  mean_raw_disc=mean_raw_disc, std_raw_disc=std_raw_disc)
>>> 0.42
```

You can also pass a precomputed similarity matrix:
```
sari = SARI(pKi=pKi, sim_matrix=precomputed_similarity_matrix)
```

### RMODI
``RMODI`` is a function and takes a distance matrix in square form, 
and a list of float, as input.

```
from rogi import RMODI
RMODI(Dx=square_dist_matrix, Y=Y)
>>> 0.42
```

The ``delta`` values used by default is ``0.625``, but can be changed with the ``delta`` argument:

```
from rogi import RMODI
RMODI(Dx=square_dist_matrix, Y=Y, delta=0.5)
>>> 0.21
```

### MODI
``MODI`` is a function and takes a distance matrix in square form, 
and a list of binary labels (`0` and `1`), as input.

```
from rogi import MODI
MODI(Dx=square_dist_matrix, Y=Y)
>>> 0.42
```

## Citation
If you make use of the ``rogi`` package in scientific publications, please cite the following article:

```
@misc{rogi,
      title={Roughness of molecular property landscapes and its impact on modellability}, 
      author={Matteo Aldeghi and David E. Graff and Nathan Frey and Joseph A. Morrone and 
              Edward O. Pyzer-Knapp and Kirk E. Jordan and Connor W. Coley},
      year={2022},
      eprint={2207.09250},
      archivePrefix={arXiv},
      primaryClass={q-bio.QM}
      }
```

If you use ``SARI``, please also cite:

```
@article{sari,
         title={SAR Index: Quantifying the Nature of Structure−Activity Relationships},
         author={Peltason, Lisa and Bajorath, J\"urgen},
         journal={J. Med. Chem.},
         publisher={American Chemical Society},
         volume={50},
         number={23},
         pages={5571--5578},
         year={2007}
         }
```

If you use ``MODI``, please also cite:

```
@article{modi,
         title={Data Set Modelability by QSAR},
         author={"Golbraikh, Alexander and Muratov, Eugene and Fourches, Denis and
                 Tropsha, Alexander"}
         journal={J. Chem. Inf. Model.},
         publisher={American Chemical Society},
         volume={54},
         number={1},
         pages={1--4},
         year={2014}
         }
```

If you use ``RMODI``, please also cite:

```
@article{rmodi,
         title={Regression Modelability Index: A New Index for Prediction of the
                Modelability of Data Sets in the Development of QSAR
                Regression Models},
         author={Luque Ruiz, Irene and G\'omez-Nieto, Miguel \'Angel},
         journal={J. Chem. Inf. Model.},
         publisher={American Chemical Society},
         volume={58},
         number={10},
         pages={2069--2084},
         year={2018}
         }
```
