#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       element_matrix.py
#       
#       Copyright 2013 Zifu Wang <z@mesh2mesh.com>
#

import numpy as np
import quadrature_class as qr
import reference_mapping as ref

def w0_w0(elt,qrule=qr.lin_tetra_2):
	""" Scalar product of node functions in a given elt, 
		return in coo_matrix format: row, col and val"""
	n_rows=elt.n_nodes
	n_cols=elt.n_nodes
	row=np.repeat(elt.nodes,n_cols)
	col=np.array(elt.nodes*n_rows)
	val=np.zeros([n_rows,n_cols],dtype=float)
	for i in xrange(qrule.n_pts):
		w0=ref.hgrad_mapping(elt,ref.w0_interpolation(elt, qrule.pts[i,:]))
		val=val+np.dot(np.transpose(w0),w0)*qrule.weights[i]*abs(elt.jacob_det)
	return row, col, np.reshape(val,-1)

def w1_w1(elt,qrule=qr.lin_tetra_2):
	""" Scalar product of edge functions in a given elt, 
		return in coo_matrix format: row, col and val"""
	n_rows=elt.n_edges
	n_cols=elt.n_edges
	row=np.repeat(elt.edges,n_cols)
	col=np.array(elt.edges*n_rows)
	val=np.zeros([n_rows,n_cols],dtype=float)
	for i in xrange(qrule.n_pts):
		w1=ref.hcurl_mapping(elt,ref.w1_interpolation(elt, qrule.pts[i,:]))
		val=val+np.dot(np.transpose(w1),w1)*qrule.weights[i]*abs(elt.jacob_det)
	w=elt.w_edges
	ww=np.dot(np.transpose(np.matrix(w)),np.matrix(w))
	val=np.multiply(val,ww)
	return row, col, np.reshape(val,-1)

def w2_w2(elt,qrule=qr.lin_tetra_2):
	""" Scalar product of facet functions in a given elt, 
		return in coo_matrix format: row, col and val"""
	n_rows=elt.n_facets
	n_cols=elt.n_facets
	row=np.repeat(elt.facets,n_cols)
	col=np.array(elt.facets*n_rows)
	val=np.zeros([n_rows,n_cols],dtype=float)
	for i in xrange(qrule.n_pts):
		w2=ref.hdiv_mapping(elt,ref.w2_interpolation(elt, qrule.pts[i,:]))
		val=val+np.dot(np.transpose(w2),w2)*qrule.weights[i]*abs(elt.jacob_det)
	w=elt.w_facets
	ww=np.dot(np.transpose(np.matrix(w)),np.matrix(w))
	val=np.multiply(val,ww)
	return row, col, np.reshape(val,-1)

def w3_w3(elt,indice,qrule=qr.lin_tetra_1):
	""" Scalar product of volume functions in a given elt, 
		return in coo_matrix format: row, col and val"""
	n_rows=1
	n_cols=1
	row=np.array([indice])
	col=np.array([indice])
	val=np.zeros([n_rows,n_cols],dtype=float)
	for i in xrange(qrule.n_pts):
		w3=ref.vol_mapping(elt,ref.w3_interpolation(elt, qrule.pts[i,:]))
		val=val+np.dot(np.transpose(w3),w3)*qrule.weights[i]*abs(elt.jacob_det)
	return row, col, np.reshape(val,-1)

	
