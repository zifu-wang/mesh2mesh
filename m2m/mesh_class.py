#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       mesh_class.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#       

import os
import sys
import bisect
import numpy as np

import m2m.io_base as io
import m2m.math_tools as mt
import m2m.simplex_class as simplex
import m2m.element_type_class as eltype

class Mesh(object):
    def __init__(self, name=None, 
        n_nodes=None, node_list=None, n_edges=None, edge_list=None, 
        n_facets=None, facet_list=None, n_elements=None, element_list=None ):
        #name
        self.name = name
        #node
        self.n_nodes = n_nodes
        if node_list==None: 
            self.node_list = []
        else:
            self.node_list = node_list
        self.map_nodes=None # irregular node labels
        #edge
        self.n_edges = n_edges
        if edge_list==None: 
            self.edge_list = []
        else:
            self.edge_list = edge_list
        #facet
        self.n_facets = n_facets
        if facet_list==None: 
            self.facet_list = []
        else:
            self.facet_list = facet_list
        #element
        self.n_elements = n_elements
        if element_list==None: 
            self.element_list = []
        else:
            self.element_list = element_list
        # kdtree
        self.tree=None

    def build_nparrays(self):
        #for each node: position, n_elts, elts
        self.np_nodes = np.array([i.pos for i in self.node_list],dtype=float)
        self.np_node2n_elts = np.array([i.n_elts for i in self.node_list],dtype=int)
        self.np_node2elts = np.empty([self.n_nodes,np.max(self.np_node2n_elts)],
                                     dtype=int)
        for inode in range(self.n_nodes):
            node=self.node_list[inode]
            self.np_node2elts[inode,:node.n_elts]=node.elts
        #for each elt: eltype, n_nodes, nodes
        self.np_elt2type = np.array([i.elt_type for i in self.element_list],dtype=int)
        self.np_elt2n_nodes = np.array([i.n_nodes for i in self.element_list],dtype=int)
        self.np_elt2nodes = np.empty([self.n_elements,eltype.MAX_N_NODES],dtype=int)
        for ielt in range(self.n_elements):
            elt=self.element_list[ielt]
            self.np_elt2nodes[ielt,:elt.n_nodes]=elt.nodes
    

def build_edges(mesh):
    """ Build edge list/connectity on the given mesh """
    print(" -- Building Edges")
    # firs loop, build keys
    keys=[]
    # copy globle function to local
    cantor_pair=mt.cantor_pair
    inv_cantor_pair=mt.inv_cantor_pair
    for elt in mesh.element_list:
        ref_elt=eltype.select_ref_elt(elt)
        for [i,j] in ref_elt.edges:
            inode=elt.nodes[i]
            jnode=elt.nodes[j]
            keys.append(cantor_pair(sorted([inode, jnode])))
    # sort keys and build edge list 
    uniq_keys=sorted(set(keys))
    mesh.n_edges=len(uniq_keys)
    for i in uniq_keys:
        mesh.edge_list.append(simplex.Edge(ident=i,nodes=inv_cantor_pair(i)))
    # second loop, update edges in elts
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        ref_elt=eltype.select_ref_elt(elt)
        # init
        elt.n_edges=ref_elt.n_edges 
        elt.edges=[]; elt.w_edges=[]
        for [i,j] in ref_elt.edges:
            inode=elt.nodes[i]
            jnode=elt.nodes[j]
            edge_key=cantor_pair(sorted([inode, jnode]))    
            edge_index=bisect.bisect_left(uniq_keys,edge_key)
            elt.edges.append(edge_index)
            mesh.edge_list[edge_index].elts.append(ielt)
            if inode<jnode: 
                elt.w_edges.append(+1)
            else:
                elt.w_edges.append(-1)
    print(io.done), 
    print(io.bcolors.INFOBLUE + '\t #: ' + str(mesh.n_edges) + io.bcolors.ENDC) 

def build_facets(mesh,only_3d=True):
    """ Build facet list/connectity on the given mesh """
    if mesh.n_facets!=None: return  # already done
    print(" -- Building Facets")
    # firs loop, build keys
    keys=[]
    # copy globle function to local
    cantor_3=mt.cantor_3
    inv_cantor_3=mt.inv_cantor_3
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if (elt.elt_type not in eltype.THREE_D_ELT_TYPES and only_3d): continue
        ref_elt=eltype.select_ref_elt(elt)
        for [i,j,k] in ref_elt.facets:
            inode=elt.nodes[i]
            jnode=elt.nodes[j]
            knode=elt.nodes[k]
            keys.append(cantor_3(sorted([inode, jnode, knode])))
    # sort keys and build facet list 
    uniq_keys=sorted(set(keys))
    mesh.n_facets=len(uniq_keys)
    for i in uniq_keys:
        mesh.facet_list.append(simplex.Facet(ident=i, nodes=inv_cantor_3(i)))
    # second loop, update facet in elts
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if (elt.elt_type not in eltype.THREE_D_ELT_TYPES and only_3d): continue
        ref_elt=eltype.select_ref_elt(elt)
        # init
        elt.n_facets=ref_elt.n_facets 
        elt.facets=[]; elt.w_facets=[]
        for [i,j,k] in ref_elt.facets:
            inode=elt.nodes[i]
            jnode=elt.nodes[j]
            knode=elt.nodes[k]
            facet_key=cantor_3(sorted([inode, jnode, knode]))
            facet_index=bisect.bisect_left(uniq_keys,facet_key)
            elt.facets.append(facet_index)
            mesh.facet_list[facet_index].on_boundary=(
                    not mesh.facet_list[facet_index].on_boundary)
            mesh.facet_list[facet_index].elts.append(ielt)
            # calculte global normal
            [xnode,ynode,znode]=sorted([inode,jnode,knode])
            xpos=mesh.node_list[xnode].pos
            ypos=mesh.node_list[ynode].pos
            zpos=mesh.node_list[znode].pos
            vxvy=ypos-xpos
            vxvz=zpos-xpos
            g_normal=np.cross(vxvy,vxvz)
            if mesh.facet_list[facet_index].center_pos==None:
                mesh.facet_list[facet_index].center_pos=(
                                    xpos+ypos+zpos)/3.
                mesh.facet_list[facet_index].normal=(g_normal/
                            np.dot(g_normal,g_normal)**0.5)
            if ([xnode,ynode,znode]==[inode,jnode,knode] or
                [xnode,ynode,znode]==[knode,inode,jnode] or
                [xnode,ynode,znode]==[jnode,knode,inode]): 
                elt.w_facets.append(+1)
            else:
                elt.w_facets.append(-1)
    print(io.done) 
    print(io.bcolors.INFOBLUE + '\t #: ' + str(mesh.n_facets) + io.bcolors.ENDC) 

def read_msh(filename,only_3d=False):
    """ Read mesh from file (*.msh).  """
    mesh=Mesh(name=os.path.splitext(filename)[0])
    #init
    line_indice=0
    import_node_mode=False
    import_node_num_mode=False
    import_element_mode=False
    import_element_num_mode=False
    
    with open(filename, 'r') as inF:
        print(io.line)
        print(io.bcolors.HEADER + 
              '    Reading mesh from file: '+ filename + io.bcolors.ENDC)
        if (only_3d): print(' [Info] 3D mode enabled: filtering all 2D elements')
        for line in inF:
            line_indice=line_indice+1
            if '$Nodes' in line:
                node_start_line=line_indice
                print('\n ++ Keyword $Nodes found in line ', line_indice)
                import_node_num_mode=True
            elif '$EndNodes' in line: 
                print(' -- Keyword $EndNodes found in line ', line_indice)
                import_node_mode=False
                #node_list check
                mesh.node_list.sort(key=lambda x: x.label)
                if (mesh.node_list[0].label!=1 or mesh.node_list[-1].label!=mesh.n_nodes):
                    print (io.bcolors.WARNING + 
                        'Node labels are not regular: thus use auto-mapping' + 
                        io.bcolors.ENDC)
                    mesh.map_nodes=[i.label for i in mesh.node_list]
            elif '$Elements' in line:
                element_start_line=line_indice
                print('\n ++ Keyword $Elements found in line ', line_indice)
                import_element_num_mode=True
            elif '$EndElements' in line: 
                print(' -- Keyword $EndElements found in line ', line_indice)
                import_element_mode=False
            elif import_node_num_mode:
                mesh.n_nodes=int(line)
                print('    - number of the nodes to be imported  = ', mesh.n_nodes)
                import_node_num_mode=False
                import_node_mode=True
            elif import_node_mode:
                line_array=[io.convert(s) for s in line.split()]
                mesh.node_list.append(simplex.Node(label=line_array[0],
                                                   x=line_array[1],
                                                   y=line_array[2],
                                                   z=line_array[3]))
            elif import_element_num_mode:
                mesh.n_elements=int(line)
                print('    - number of the elements to be imported  = ', mesh.n_elements)
                import_element_num_mode=False
                import_element_mode=True
            elif import_element_mode:
                line_array=[io.convert(s) for s in line.split()]
                if line_array[1] not in eltype.KNOWN_ELT_TYPES:
                    print(io.bcolors.WARNING + 
                    '/!\ Skipping unknown element type' + io.bcolors.ENDC)
                    print ("\tline:",line_indice,
                            "\telt label:",line_array[0],
                            "\telt type:",line_array[1]) 
                    continue
                elif only_3d and line_array[1] not in eltype.THREE_D_ELT_TYPES:
                    continue
                n_tags=line_array[2]
                mesh.element_list.append(simplex.Element(label=line_array[0],
                                                         elt_type=line_array[1],
                                                         tags=line_array[3:3+n_tags],
                                                         l_nodes=line_array[3+n_tags:],
                                                         map_nodes=mesh.map_nodes,
                                                         n_nodes=len(line_array[3+n_tags:]
                                                         )))
            
        if (only_3d): mesh.n_elements=len(mesh.element_list)
    return mesh
    
