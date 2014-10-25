#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       quadrature_class.py
#       
#       Copyright 2013 Zifu Wang <zifu.wang@icloud.com>
#

import numpy as np


class Quadrature(object):
        def __init__(self, n_pts=None, pts=None, weights=None):
		self.n_pts=n_pts
		self.pts=np.array(pts,dtype=float)
		self.weights=np.array(weights,dtype=float)

# linear tetrahedron
lin_tetra_1=Quadrature(	n_pts=1, 
			pts=[[0.25,0.25,0.25]], 
			weights=[1./6.]
			)

lin_tetra_2=Quadrature(	n_pts=4, 
			pts=[[0.138196601125011,0.138196601125011,0.138196601125011],
			     [0.138196601125011,0.138196601125011,0.585410196624968],
			     [0.138196601125011,0.585410196624968,0.138196601125011],
			     [0.585410196624968,0.138196601125011,0.138196601125011]],
			weights=[1./24.,1./24.,1./24.,1./24.]
			)
