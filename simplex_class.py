#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       simplex_class.py
#       
#       Copyright 2013 Zifu Wang <z@mesh2mesh.com>
#       

import numpy as np
import bisect

class Node(object):
	def __init__(self, label=None, x=None, y=None, z=None):
		self.label = label	 #/!\ very often different from index
		self.x = x
		self.y = y
		self.z = z
		self.pos = np.array([x, y, z],dtype=float)

class Edge(object):
	def __init__(self, ident=None, nodes=None):
		self.ident = ident	#Cantor pairing calculated from nodes
		self.nodes = nodes	#node index
		self.elts=[]		#elt index

class Facet(object):
	def __init__(self,  ident=None, nodes=None):
		self.ident = ident	#Cantor pairing calculated from nodes
		self.nodes = nodes      #node index
		self.elts=[]		#elt index
		
		self.center_pos = None 	
		self.normal = None	
		self.on_boundary=False	

class Element(object):
	def __init__(self, label=None, elt_type=None,
			n_nodes=None, nodes=None, l_nodes=None, map_nodes=None,
			n_edges=None, edges=None, w_edges=None,  
			n_facets=None, facets=None, w_facets=None,
			tags=None ):
				
		self.label = label	 #/!\ very often different from index
		self.elt_type = elt_type
		
		self.n_nodes = n_nodes
		if nodes != None: 
			self.nodes = nodes	#node index
		elif l_nodes!=None and map_nodes==None: 	#node labels -> node index
			self.nodes = [x-1 for x in l_nodes]
		elif l_nodes!=None and map_nodes!=None: 	#irrugular node labels 
			self.nodes = [bisect.bisect_left(map_nodes,x) for x in l_nodes]
		self.n_edges = n_edges
		self.edges = edges	#edge index
		self.w_edges = w_edges
		self.n_facets = n_facets
		self.facets = facets	#facet index
		self.w_facets = w_facets
		
		self.tags = tags

		self.nodes_pos = None
		self.center_pos = None
		self.jacob = None
		self.jacob_det = None
		self.jacob_inv = None

#class Group(object):
#	def __init__(self ):
#		self.num = 0
#		self.mat_list = []
#		self.elt_list = []



