#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       mesh_tree.py
#       
#       Copyright 2013 Zifu Wang <zifu.wang@icloud.com>
#       

import numpy as np

import io_base as io
import math_tools as mt
import mesh_class as mc
import simplex_class as simplex
import reference_mapping as ref
import element_type_class as eltype

from random import choice
from scipy.spatial import cKDTree

from fortran_walk import walking_search

def pre_process(mesh,want_edges=False,want_facets=False,want_tree=False):
	""" Renumber nodes, update jacob, build edge/facet list and KDtree """
	print io.line
	print (io.bcolors.HEADER + '    Pre-processing mesh: ' + mesh.name + '\n' 
							+ io.bcolors.ENDC)
	renumber_nodes(mesh)
	if want_edges: mc.build_edges(mesh)
	if want_facets or want_tree: mc.build_facets(mesh)
	if want_tree: build_tree(mesh)

def renumber_nodes(mesh,only_3d=True):
	""" Renumber nodes in all elements and initialize jacob matrix"""
	print ' -- Updating Nodes & Jacobian matrix'
        for ielt in xrange(mesh.n_elements):
                elt=mesh.element_list[ielt]
		if (elt.elt_type !=4): continue
		#linear tetra
		ref.update_jacob(elt, mesh.node_list)
		if elt.jacob_det > 0: continue
		(elt.nodes[1],elt.nodes[2])=(elt.nodes[2],elt.nodes[1])
		ref.update_jacob(elt, mesh.node_list, forced=True)
	print io.done

def build_tree(mesh,leafs=12,only_3d=True):
	""" Build cKDTree on the given mesh """
	if mesh.tree!=None: return # already done
	print ' -- Building cKDTree'
        bary_list=np.empty([mesh.n_elements,3],dtype=float) 
        for ielt in xrange(mesh.n_elements):
                elt=mesh.element_list[ielt]
		if (elt.elt_type not in eltype.THREE_D_ELT_TYPES and only_3d): 
			bary_list[ielt,:]=None
		else:
                	bary_list[ielt,:]=elt.center_pos[:]
        mesh.tree=cKDTree(bary_list,leafsize=leafs)
	print io.done

def fort_find_containing_elt(mesh,pts):
	""" Find the pt-containing elt using cKDTree + walking """
	n_pts=pts.shape[0]
	print io.bcolors.INFOBLUE + '\t+ locating: ', n_pts, ' points', 
	print ' on mesh: ' + mesh.name + io.bcolors.ENDC 
	print '\t|-- Querying cKDTree'
	(d,inits)=mesh.tree.query(pts,k=1,eps=0.5)
	print '\t|' + io.done
	print '\t|-- Preparing fortran interface'+ io.bcolors.ENDC
	inits=np.array(inits,dtype=int)+1
	n_elts=mesh.n_elements
	n_facets=mesh.n_facets
	elt_max_n_facets=max([i.n_facets for i in mesh.element_list])
	elt_n_facets=np.array([i.n_facets for i in mesh.element_list],dtype=int)
	elt_facets=np.empty([n_elts,elt_max_n_facets],dtype=int)
	elt_w_facets=np.empty([n_elts,elt_max_n_facets],dtype=int)
	elt_outgoings=np.empty([n_elts,elt_max_n_facets],dtype=int)
        for ielt in xrange(mesh.n_elements):
                elt=mesh.element_list[ielt]
		ref_elt=eltype.select_ref_elt(elt)
		if elt.n_facets==None: continue
		elt_facets[ielt,:elt.n_facets]=np.array(elt.facets,dtype=int)+1	
		elt_w_facets[ielt,:elt.n_facets]=elt.w_facets
		elt_outgoings[ielt,:elt.n_facets]=ref_elt.outgoing
	facet_pos=np.empty([n_facets,3],dtype=float)
	facet_normals=np.empty([n_facets,3],dtype=float)
	facet_elts=np.empty([n_facets,2],dtype=int)
        for ifacet in xrange(mesh.n_facets):
		facet=mesh.facet_list[ifacet]
		facet_pos[ifacet,:]=facet.center_pos
		facet_normals[ifacet,:]=facet.normal
		facet_elts[ifacet,:]=-1
		facet_elts[ifacet,:len(facet.elts)]=np.array(facet.elts,dtype=int)+1
	print '\t|' + io.done
	print '\t|-- Walking thru elements (fortran_walk.so)'
	res=walking_search(elt_n_facets,elt_facets,elt_w_facets,elt_outgoings,
				facet_pos,facet_normals,facet_elts,pts,inits,99,
				n_elts,n_facets,elt_max_n_facets,n_pts)
	res=np.array(res,dtype=int)-1
	print '\t~' + io.done
	return res
	
