#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random
import cProfile

import numpy as np
import mesh_class 
import mesh_tree
import reference_mapping as ref

import field_projection

filename="test8.msh"
mesh=mesh_class.read_msh(filename)

mesh_tree.pre_process(mesh,want_edges=True,want_facets=True,want_tree=True)

#print mesh_tree.is_in_elt(mesh.element_list[2551],[0.1875,  0.8125,  0.4376])
#for i in range(1000):
	#print mesh_tree.fort_find_containing_elt(mesh,[ 0.1875,  0.8125,  0.4376])

	#print mesh_tree.fort_walking_search_containing_elt(mesh,[ 0.1875,  0.8125,  0.4376],init=100)
#pts=np.empty([100000,3],dtype=float)
#for i in xrange(100000):
#	pt=np.array([random.random(),random.random(),random.random()])
#	pts[i,:]=pt
print mesh_tree.fort_find_containing_elt(mesh,np.array([[0.45,0.35,0.25],[0.5,0.6,0.7]]))

#print node info    
#for inode in range(mesh.n_nodes):
#	node=mesh.node_list[inode]
#	print inode, node.label 
#	print node.pos 

#print elt info    
#for ielt in range(mesh.n_elements):
#	elt=mesh.element_list[ielt]
#	print ielt
#	print elt.elt_type,elt.n_nodes,elt.nodes 

#field=np.empty(mesh.n_elements,dtype=float)
#for ielt in range(mesh.n_elements):
#	elt=mesh.element_list[ielt]
#	field[ielt]=elt.center_pos[2]**2
##	print ielt, elt.label, elt.center_pos, elt.jacob_det
##	print elt.nodes 
##	print elt.edges,elt.w_edges 
##	print elt.facets,elt.w_facets 
#np.savetxt("sfield.csv", field)
#for iedge in range(mesh.n_edges):
#        edge=mesh.edge_list[iedge]
#	print iedge, 
#	print edge.nodes, edge.elts
#
#for ifac in range(mesh.n_facets):
#        facet=mesh.facet_list[ifac]
#	print ifac, 
#	print facet.nodes, facet.normal, facet.on_boundary, facet.elts

#field_projection.main(mesh)
