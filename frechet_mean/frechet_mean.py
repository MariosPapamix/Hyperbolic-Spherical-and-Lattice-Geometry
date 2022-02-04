import os
import sys
import warnings
import numpy as np

sys.path.append(os.path.dirname(os.getcwd()))
warnings.filterwarnings('ignore')

import geomstats.backend as gs
from geomstats.learning.frechet_mean import FrechetMean
from geomstats.geometry.hypersphere import Hypersphere
from geomstats.geometry.poincare_ball import PoincareBall

def calc_frechet_mean(data, d, geom, maxit):
    d = int(d) # make sure int
    # initialise geometry
    if geom == "poincare_disk":
        space = PoincareBall(dim=d)
    elif geom == "sphere":
        space = Hypersphere(dim=d)
    # calculate frechet mean   
    mean = FrechetMean(metric=space.metric)
    mean.max_iter = maxit
    mean.fit(data)
    est = mean.estimate_
    if geom == "sphere":
        est = est / sum(est) # noramlise to be on unit sphere
    return est
