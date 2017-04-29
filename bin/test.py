#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random
import cProfile

import numpy as np
import m2m.mesh_class as mc
import m2m.mesh_tree as mt

#import m2m.field_projection as fp

filename="test32.msh"
mesh=mc.read_msh(filename,only_3d=True)
mt.pre_process(mesh,want_edges=True,want_facets=True,want_tree=True)
#cProfile.run('mt.pre_process(mesh,want_edges=True,want_facets=True,want_tree=True)')

#np.savetxt("node2pos.csv", mesh.np_nodes, fmt='%e')
#np.savetxt("elt2nodes.csv", mesh.np_elt2nodes, fmt='%i')
#np.savetxt("elt2type.csv", mesh.np_elt2type, fmt='%i')
#np.savetxt("facet2cpos.csv", mesh.np_facet2cpos, fmt='%e')
#np.savetxt("facet2og.csv", mesh.np_facet2og, fmt='%e')
#np.savetxt("elt2facets.csv", mesh.np_elt2facets, fmt='%i')
#np.savetxt("facet2elts.csv", mesh.np_facet2elts, fmt='%i')
#np.savetxt("elt2w_facets.csv", mesh.np_elt2w_facets, fmt='%i')

print(mt.fort_find_containing_elt(mesh,np.array([[0.45,0.35,0.25],[0.5,0.6,0.7]])))

