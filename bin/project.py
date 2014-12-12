#!/usr/bin/env python
# -*- coding: utf-8 -*-
import profile
import numpy as np

import mesh_class 
import mesh_tree 
import field_projection

smesh_filename="test8.msh"
sfield_filename="test8_elt_x2.csv"
tmesh_filename="test.msh"
target_simplex="edge"

smesh=mesh_class.read_msh(smesh_filename)
mesh_tree.pre_process(smesh,want_edges=True,want_facets=True,want_tree=True)
sfield=np.loadtxt(sfield_filename)

tmesh=mesh_class.read_msh(tmesh_filename)
mesh_tree.pre_process(tmesh,want_edges=True,want_facets=True,want_tree=False)

field_projection.main(smesh,sfield,tmesh,tdof=target_simplex)
