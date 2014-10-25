#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       element_type_class.py
#       
#       Copyright 2013 Zifu Wang <zifu.wang@icloud.com>
#

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
                self.edges = edges      
                self.n_facets = n_facets
                self.facets = facets    
                self.outgoing = outgoing 

KNOWN_ELT_TYPES=[
		1,	#linear edge
                2,      #linear triangle
                4       #linear tetrahedron
                ]

THREE_D_ELT_TYPES=[
 		  4     #linear tetrahedron
                  ]

def select_ref_elt(elt):
	if elt.elt_type==1: return lin_edge
	elif elt.elt_type==2: return lin_trian
	elif elt.elt_type==4: return lin_tetra

lin_edge=Elt_type( name="Linear edge",
		    elt_type=1,
		    n_nodes=2,
		    n_edges=1,
		    edges=[[0,1]],
                    n_facets=0,
		    facets=[],
		    outgoing=None
                  )

lin_trian=Elt_type( name="Linear triangle",
		    elt_type=2,
		    n_nodes=3,
		    n_edges=3,
		    edges=[[0,1],[0,2],[1,2]],
                    n_facets=1,
		    facets=[[0,1,2]],
		    outgoing=None
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


