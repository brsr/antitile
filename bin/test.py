# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:53:50 2017

@author: Bstone
"""

import numpy as np
import antitile

sq = antitile.xmath.normalize(np.array([[1,2,1],
                               [-1,2,1],
                               [-1,-1,1],
                               [1,-2,1]]))
bkdn = antitile.breakdown.Breakdown(16,0,shape=4)
xy = bkdn.coord[:, np.newaxis]
x1 = antitile.projection.square_slerp(xy, sq)
x2 = antitile.projection.square_slerp(xy[:,::-1], sq[[0,3,2,1]])