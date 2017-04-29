#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       mesh_tree.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#       

import numpy as np

import m2m.io_base as io
import m2m.reference_mapping as ref
import m2m.element_type_class as eltype

import m2m.lib.fortran_topo as ftopo
import m2m.lib.fortran_walk as fwalk

from scipy.spatial import cKDTree


def pre_process(mesh,want_edges=False,want_facets=False,want_tree=False):
    """ Calculate jacob matrix, build edge/facet list and KDtree """
    print(io.line)
    print(io.bcolors.HEADER + '    Pre-processing mesh: ' + mesh.name + '\n' 
                            + io.bcolors.ENDC)
    # build node2elts list
    build_node2elts(mesh)
    # build numpy arrays: faster slicing & call fortran subroutines
    mesh.build_nparrays()
    # update node position in each elt and calculate Jacobian stuff
    calculate_eltpos_jacob(mesh)
    if want_edges: build_edges(mesh) 
    if want_facets: build_facets(mesh)
    if want_tree: build_tree(mesh)

def build_node2elts(mesh):
    """ Build element list for each node """
    print(' -- Building node2elts list')
    for ielt in range(mesh.n_elements):
        for node in mesh.element_list[ielt].nodes:
            mesh.node_list[node].elts.append(ielt)
            mesh.node_list[node].n_elts=mesh.node_list[node].n_elts+1
    print(io.done)

def calculate_eltpos_jacob(mesh,renumber=True):
    """ Calculate jacob matrix """
    print(' -- Updating element position & Jacobian stuff')
    #if renumber: print(' [Info] Renumber mode enabled: renumber nodes if negative J-det')
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        elt.nodes_pos=mesh.np_nodes[mesh.np_elt2nodes[ielt,:mesh.np_elt2n_nodes[ielt]]]
        (elt.center_pos, elt.jacob, elt.jacob_det, elt.jacob_inv, elt.det_sign) = \
                                ref.update_jacob(elt.elt_type, elt.nodes_pos)
        if renumber and elt.det_sign<0:
            ref.renumber_elt(elt.elt_type,mesh,ielt)
            elt.nodes_pos = \
                    mesh.np_nodes[mesh.np_elt2nodes[ielt,:mesh.np_elt2n_nodes[ielt]]]
            (elt.center_pos, elt.jacob, elt.jacob_det, elt.jacob_inv, elt.det_sign) = \
                                ref.update_jacob(elt.elt_type, elt.nodes_pos)
    print(io.done)

def build_edges(mesh):
    """ Build edge list/connectity on the given mesh """
    # only work for tetra
    print(" -- Building Edges",end="")
    (mesh.n_edges,mesh.elt2max_n_edges,mesh.np_elt2n_edges,mesh.np_node2n_edges) = \
        ftopo.count_edges(elt2n_nodes=mesh.np_elt2n_nodes,elt2nodes=mesh.np_elt2nodes,
                         node2n_elts=mesh.np_node2n_elts,node2elts=mesh.np_node2elts)
    print("\t",mesh.n_edges," edges found")
    (mesh.np_edge2nodes,mesh.np_elt2edges,mesh.np_elt2w_edges) = \
        ftopo.build_edges(elt2n_nodes=mesh.np_elt2n_nodes,elt2nodes=mesh.np_elt2nodes,
                          node2n_elts=mesh.np_node2n_elts,node2elts=mesh.np_node2elts,
                          n_edges=mesh.n_edges,elt2max_n_edges=mesh.elt2max_n_edges,
                          elt2n_edges=mesh.np_elt2n_edges)
    print(io.done)

def build_facets(mesh):
    """ Build facet list/connectity on the given mesh """
    # only work for tetra
    print(" -- Building Facets",end="")
    (mesh.n_facets,mesh.elt2max_n_facets,mesh.np_elt2n_facets,mesh.np_node2n_facets) = \
        ftopo.count_facets(elt2n_nodes=mesh.np_elt2n_nodes,elt2nodes=mesh.np_elt2nodes,
                         node2n_elts=mesh.np_node2n_elts,node2elts=mesh.np_node2elts)
    print("\t",mesh.n_facets," facets found")
    (mesh.np_facet2nodes,mesh.np_facet2elts,mesh.np_facet2cpos,mesh.np_facet2normal,
        mesh.np_elt2facets,mesh.np_elt2w_facets) = \
        ftopo.build_facets(elt2n_nodes=mesh.np_elt2n_nodes,elt2nodes=mesh.np_elt2nodes,
                          node2n_elts=mesh.np_node2n_elts,node2elts=mesh.np_node2elts,
                          n_facets=mesh.n_facets,elt2max_n_facets=mesh.elt2max_n_facets,
                          elt2n_facets=mesh.np_elt2n_facets,nodes=mesh.np_nodes)
    print(io.done)

def build_tree(mesh,leafs=12):
    """ Build cKDTree on the given mesh """
    if mesh.tree!=None: return # already done
    print(' -- Building cKDTree')
    bary_list=np.array([i.center_pos for i in mesh.element_list],dtype=float) 
    mesh.tree=cKDTree(bary_list,leafsize=leafs)
    print(io.done)

def fort_find_containing_elt(mesh,pts):
    """ Find the pt-containing elt using cKDTree + walking """
    print('\t+ locating: ', pts.shape[0], ' points on mesh: ', mesh.name)
    print('\t|-- Querying cKDTree')
    (d,inits)=mesh.tree.query(pts,k=1,eps=0.5)
    print('\t|' + io.done)
    print('\t|-- Preparing fortran interface'+ io.bcolors.ENDC)
    inits=np.array(inits,dtype=int)
    mesh.np_elt2facet_og=np.empty([mesh.n_elements,mesh.elt2max_n_facets],dtype=int)
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        ref_elt=eltype.select_ref_elt(elt)
        mesh.np_elt2facet_og[ielt,:elt.n_facets]=ref_elt.outgoing
    print('\t|' + io.done)
    print('\t|-- Walking thru elements (multi-processing fortran_walk.so)')
    res=fwalk.walking_search(elt_n_facets=mesh.np_elt2n_facets,
                             elt_facets=mesh.np_elt2facets,
                             elt_w_facets=mesh.np_elt2w_facets,
                             elt_outgoings=mesh.np_elt2facet_og,
                             facet_pos=mesh.np_facet2cpos,
                             facet_normals=mesh.np_facet2normal,
                             facet_elts=mesh.np_facet2elts,
                             pts=pts,inits=inits,maxtry=99)
    print('\t~' + io.done)
    return res
    
