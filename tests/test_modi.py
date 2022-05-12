#!/usr/bin/env python
from rogi import MODI, RMODI
import math
import numpy as np
from scipy.spatial.distance import squareform


def test_modi():
    Y = [0., 1., 0., 1., 0.]
    d = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]) / 10.
    Dx = squareform(d)
    score = MODI(Dx=Dx, Y=Y)
    assert math.isclose(score, 0.3333333333333333)


def test_rmodi():
    Y = [1., 3., 2., 1., 4.]
    d = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]) / 10.
    Dx = squareform(d)
    score = RMODI(Dx=Dx, Y=Y)
    assert math.isclose(score, 0.2)
