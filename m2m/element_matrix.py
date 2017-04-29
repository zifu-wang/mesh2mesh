#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       element_matrix.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#

import numpy as np
import m2m.quadrature_class as qr
import m2m.reference_mapping as ref

def w0_w0(mesh,ielt,qrule=qr.lin_tetra_2):
    """ Scalar product of node functions in a given elt, 
        return in coo_matrix format: row, col and val"""
    elt=mesh.element_list[ielt]
    n_rows=mesh.np_elt2n_nodes[ielt]
    n_cols=mesh.np_elt2n_nodes[ielt]
    row=np.repeat(mesh.np_elt2nodes[ielt,:],n_cols)
    col=np.tile(mesh.np_elt2nodes[ielt,:],n_rows)
    val=np.zeros([n_rows,n_cols],dtype=float)
    for i in range(qrule.n_pts):
        w0=ref.hgrad_mapping(elt,ref.w0_interpolation(mesh,ielt,qrule.pts[i,:]))
        val=val+np.dot(np.transpose(w0),w0)*qrule.weights[i]*abs(elt.jacob_det)
    return row, col, np.reshape(val,-1)

def w1_w1(mesh,ielt,qrule=qr.lin_tetra_2):
    """ Scalar product of edge functions in a given elt, 
        return in coo_matrix format: row, col and val"""
    elt=mesh.element_list[ielt]
    n_rows=mesh.np_elt2n_edges[ielt]
    n_cols=mesh.np_elt2n_edges[ielt]
    row=np.repeat(mesh.np_elt2edges[ielt,:],n_cols)
    col=np.tile(mesh.np_elt2edges[ielt,:],n_rows)
    val=np.zeros([n_rows,n_cols],dtype=float)
    for i in range(qrule.n_pts):
        w1=ref.hcurl_mapping(elt,ref.w1_interpolation(mesh,ielt,qrule.pts[i,:]))
        val=val+np.dot(np.transpose(w1),w1)*qrule.weights[i]*abs(elt.jacob_det)
    w=mesh.np_elt2w_edges[ielt,:]
    ww=np.dot(np.transpose(np.matrix(w)),np.matrix(w))
    val=np.multiply(val,ww)
    return row, col, np.reshape(val,-1)

def w2_w2(mesh,ielt,qrule=qr.lin_tetra_2):
    """ Scalar product of facet functions in a given elt, 
        return in coo_matrix format: row, col and val"""
    elt=mesh.element_list[ielt]
    n_rows=mesh.np_elt2n_facets[ielt]
    n_cols=mesh.np_elt2n_facets[ielt]
    row=np.repeat(mesh.np_elt2facets[ielt,:],n_cols)
    col=np.tile(mesh.np_elt2facets[ielt,:],n_rows)
    val=np.zeros([n_rows,n_cols],dtype=float)
    for i in range(qrule.n_pts):
        w2=ref.hdiv_mapping(elt,ref.w2_interpolation(mesh,ielt,qrule.pts[i,:]))
        val=val+np.dot(np.transpose(w2),w2)*qrule.weights[i]*abs(elt.jacob_det)
    w=mesh.np_elt2w_facets[ielt,:]
    ww=np.dot(np.transpose(np.matrix(w)),np.matrix(w))
    val=np.multiply(val,ww)
    return row, col, np.reshape(val,-1)

def w3_w3(mesh,ielt,qrule=qr.lin_tetra_1):
    """ Scalar product of volume functions in a given elt, 
        return in coo_matrix format: row, col and val"""
    elt=mesh.element_list[ielt]
    n_rows=1
    n_cols=1
    row=np.array([ielt])
    col=np.array([ielt])
    val=np.zeros([n_rows,n_cols],dtype=float)
    for i in range(qrule.n_pts):
        w3=ref.vol_mapping(elt,ref.w3_interpolation(mesh,ielt,qrule.pts[i,:]))
        val=val+np.dot(np.transpose(w3),w3)*qrule.weights[i]*abs(elt.jacob_det)
    return row, col, np.reshape(val,-1)

    
