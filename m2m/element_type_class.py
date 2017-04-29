#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       element_type_class.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#

import numpy as np

class Elt_type(object):
     def __init__(self, name=None, elt_type=None,
                  n_nodes=None, 
                  n_edges=None, edges=None, 
                  n_facets=None, facets=None, outgoing=None
                 ):
        
        self.name = name
        self.elt_type = elt_type
        
        self.n_nodes = n_nodes
        self.n_edges = n_edges
        self.n_facets = n_facets
        if edges == None: self.edges = edges
        else: self.edges = np.array(edges, dtype=int)      
        if facets == None: self.facets = facets
        else: self.facets = np.array(facets, dtype=int)    
        if outgoing == None: self.outgoing = outgoing
        else: self.outgoing = np.array(outgoing, dtype=int)

lin_edge=Elt_type( name="Linear edge",
                   elt_type=1,
                   n_nodes=2,
                   n_edges=1,
                   edges=[[0,1]],
                   n_facets=0,
                   facets=[],
                   outgoing=[]
                  )

lin_trian=Elt_type( name="Linear triangle",
                    elt_type=2,
                    n_nodes=3,
                    n_edges=3,
                    edges=[[0,1],[0,2],[1,2]],
                    n_facets=1,
                    facets=[[0,1,2]],
                    outgoing=[]
                  )

lin_tetra=Elt_type( name="Linear tetrahedron",
                    elt_type=4,
                    n_nodes=4,
                    n_edges=6,
                    edges=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]],
                    n_facets=4,
                    facets=[[0,1,2],[0,1,3],[0,2,3],[1,2,3]],
                    outgoing=[-1,+1,-1,+1]
                  )

def select_ref_elt(elt):
    if elt.elt_type==1: return lin_edge
    elif elt.elt_type==2: return lin_trian
    elif elt.elt_type==4: return lin_tetra
    else: return None
    
def select_ref_elt_by_type(elt_type):
    if elt_type==1: return lin_edge
    elif elt_type==2: return lin_trian
    elif elt_type==4: return lin_tetra
    else: return None

KNOWN_ELT_TYPES=[
        1,      #linear edge
        2,      #linear triangle
        4       #linear tetrahedron
                ]

THREE_D_ELT_TYPES=[
        4       #linear tetrahedron
                ]

MAX_N_NODES=max([select_ref_elt_by_type(i).n_nodes for i in KNOWN_ELT_TYPES])

MAX_N_EDGES=max([select_ref_elt_by_type(i).n_edges for i in KNOWN_ELT_TYPES])

MAX_N_FACETS=max([select_ref_elt_by_type(i).n_facets for i in KNOWN_ELT_TYPES])


