#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       reference_mapping.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#

import numpy as np
import m2m.lib.fortran_jmap as jmap

def update_jacob(elt_type,nodes_pos):
    """ Update of the Jacobian stuff """
    if elt_type == 4:   #tetrahedron
        return jmap.tetra_jacob(nodes_pos)

def renumber_elt(elt_type,mesh,ielt):
    """ renumber nodes in an elt in order to change the sign of its J-det """
    if elt_type == 4:   #tetrahedron
        mesh.np_elt2nodes[ielt,[1,2]]=mesh.np_elt2nodes[ielt,[2,1]]

def get_real_pos(elt,ref_pt):
    """ Return the coordinate position of a reference point """
    #tetrahedron
    pos=np.array(elt.nodes_pos[0,:]
            +ref_pt[0]*elt.jacob[0,:]
            +ref_pt[1]*elt.jacob[1,:]
            +ref_pt[2]*elt.jacob[2,:],
             dtype=float)
    return pos

def get_ref_pos(elt,pt):
    """ Return the reference position of a real point """
    #tetrahedron
    v=np.arrary(pt-elt.nodes_pos[0,:],dtype=float)
    ref_pos=np.dot(v,elt.jacob_inv)
    return ref_pos

def w0_interpolation(mesh,ielt,ref_pt):
    """ Return the values of REF node functions on a reference point """
    values=np.empty([1,mesh.np_elt2n_nodes[ielt]],dtype=float)
    [x,y,z]=ref_pt
    #linear tetrahedron
    values[0,0]=1-x-y-z
    values[0,1]=x
    values[0,2]=y
    values[0,3]=z
    return values

def w0grad_interpolation(mesh,ielt,ref_pt):
    """ Return the values of the gradients of
        REF node functions on a reference point """
    #linear tetrahedron
    values=np.array([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]],dtype=float)
    return values

def w1_interpolation(mesh,ielt,ref_pt):
    """ Return the values of REF edge functions on a reference point """
    values=np.empty([3,mesh.np_elt2n_edges[ielt]],dtype=float)
    [x,y,z]=ref_pt
    #linear tetrahedron
    values[:,0]=[1-y-z,x,x]
    values[:,1]=[y,1-x-z,y]
    values[:,2]=[z,z,1-x-y]
    values[:,3]=[-y,x,0]
    values[:,4]=[-z,0,x]
    values[:,5]=[0,-z,y]
    return values

def w1curl_interpolation(mesh,ielt,ref_pt):
    """ Return the values of the curls of
        REF edge functions on a reference point """
    #linear tetrahedron
    values=np.array([[0,2,-2,0,0,2],[-2,0,2,0,-2,0],[2,-2,0,2,0,0]],dtype=float)
    return values

def w2_interpolation(mesh,ielt,ref_pt):
    """ Return the values of REF facet functions on a reference point """
    values=np.empty([3,mesh.np_elt2n_facets[ielt]],dtype=float)
    [x,y,z]=ref_pt
    #linear tetrahedron
    values[:,0]=[-x,-y,1-z]
    values[:,1]=[x,-1+y,z]
    values[:,2]=[1-x,-y,-z]
    values[:,3]=[x,y,z]
    values=2.*values
    return values

def w2div_interpolation(mesh,ielt,ref_pt):
    """ Return the values of the divergences of
        REF facet functions on a reference point """
    #linear tetrahedron
    values=np.array([[-6,6,-6,6]],dtype=float)
    return values

def w3_interpolation(mesh,ielt,ref_pt):
    """ Return the values of REF volume functions on a reference point """
    #linear tetrahedron
    values=np.array([[6]],dtype=float)
    return values

def hgrad_mapping(elt,ref_values):
    """ Map vectors from ref_element to real_element """
    real_values=np.copy(ref_values)
    return real_values

def hcurl_mapping(elt,ref_values):
    """ Map vectors from ref_element to real_element """
    real_values=np.dot(elt.jacob_inv, ref_values)
    return real_values

def hdiv_mapping(elt,ref_values):
    """ Map vectors from ref_element to real_element """
    real_values=np.dot(np.transpose(elt.jacob), ref_values)/elt.jacob_det
    return real_values

def vol_mapping(elt,ref_values):
    """ Map vectors from ref_element to real_element """
    real_values=ref_values/abs(elt.jacob_det)
    return real_values

