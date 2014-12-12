#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       field_projection.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#

import sys
import numpy as np

import m2m.io_base as io
import m2m.mesh_tree as mt
import m2m.quadrature_class as qr
import m2m.reference_mapping as ref
import m2m.element_matrix as eltm

import scipy.sparse
import scipy.sparse.linalg

def main(smesh,sfield,mesh,tdof=None):
    """ Mesh-to-Mesh field projection """
    """ Note: meshes should be pre_processed before calling """
    print(io.line)
    print(io.bcolors.HEADER,'    Projecting field on: ', tdof,' dofs\n', io.bcolors.ENDC)
    print(' -- Preparing memory')
    if tdof=="node":
        cmat_n_fill=mesh.n_elements*4**2
        n_dofs=mesh.n_nodes
    elif tdof=="edge":
        cmat_n_fill=mesh.n_elements*6**2
        n_dofs=mesh.n_edges
    elif tdof=="facet":
        cmat_n_fill=mesh.n_elements*4**2
        n_dofs=mesh.n_facets
    elif tdof=="element":
        cmat_n_fill=mesh.n_elements
        n_dofs=mesh.n_elements
    else:
        sys.exit(io.bcolors.FAIL + 'Wrong target dof ' + io.bcolors.ENDC)
    cmat_row=np.empty(cmat_n_fill,dtype=int)
    cmat_col=np.empty(cmat_n_fill,dtype=int)
    cmat_val=np.empty(cmat_n_fill,dtype=float)
    fvec=np.zeros(n_dofs,dtype=float)
    print(io.done)

    print(' -- Locating sample points')
    # first loop, count number
    n_smpl_pts=0
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if elt.elt_type!=4: continue
        q_rule=qr.lin_tetra_2
        n_smpl_pts=n_smpl_pts+q_rule.n_pts
    pts=np.empty([n_smpl_pts,3],dtype=float)
    # second loop, stock pts
    i_pt=0
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if elt.elt_type!=4: continue
        for iq in range(q_rule.n_pts):
            pts[i_pt,:]=(ref.get_real_pos(elt,q_rule.pts[iq,:]))
            i_pt=i_pt+1
    # now call the searching routine
    pts2elts=mt.fort_find_containing_elt(smesh,pts)

    print(' -- Assemblying matrix')
    vol=0.
    i_pt=0
    cmat_index=0
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if elt.elt_type!=4: continue

        ## Cmat
        if tdof=="node":
            cmat_grow=4**2
            (row,col,val)=eltm.w0_w0(mesh,ielt)
        elif tdof=="edge":
            cmat_grow=6**2
            (row,col,val)=eltm.w1_w1(mesh,ielt)
        elif tdof=="facet":
            cmat_grow=4**2
            (row,col,val)=eltm.w2_w2(mesh,ielt)
        elif tdof=="element":
            cmat_grow=1
            (row,col,val)=eltm.w3_w3(mesh,ielt)

        #assembly
        cmat_row[cmat_index:cmat_index+cmat_grow]=row
        cmat_col[cmat_index:cmat_index+cmat_grow]=col
        cmat_val[cmat_index:cmat_index+cmat_grow]=val
        cmat_index=cmat_index+cmat_grow

        #node biort
#        for i in range(elt.n_nodes):
#            inode=elt.nodes[i]
#            cmat[inode,0]=cmat[inode,0]+abs(elt.jacob_det)
        #elt
#        cmat[ielt,0]=6./abs(elt.jacob_det)
        #q_rule=quadrature_class.lin_tetra_2
        #for iq=range(q_rule.n_pts)        
        #    w_elt=6/abs(elt.jacob_det)
        #    cmat[ielt]=cmat[ielt]+
        #        w_elt**2*q_rule.weights[iq]*abs(elt.jacob_det)
        
        ## Fvec
        q_rule=qr.lin_tetra_2
        for iq in range(q_rule.n_pts):
            # Quadrature pt
            real_pos=ref.get_real_pos(elt,q_rule.pts[iq,:])
            # locate on source mesh
            iselt= pts2elts[i_pt]
            i_pt=i_pt+1
            if iselt==-1: continue
            # field
#            field=real_pos[0]**2    #node,volume test
#            field=np.array([real_pos[0]**2,0,0],dtype=float)

            field=[0,sfield[iselt],0]    # todo: vector field & interpolation
#            field=sfield[iselt,-3:]    # todo: vector field & interpolation
#            field=sfield[iselt]    # todo: vector field & interpolation

            if tdof=="node":
                w0=ref.hgrad_mapping(elt,
                    ref.w0_interpolation(mesh,ielt,q_rule.pts[iq,:]))
                for i in range(mesh.np_elt2n_nodes[ielt]):
                    inode=mesh.np_elt2nodes[ielt,i]
                    fvec[inode]=(fvec[inode]+     
                            w0[0,i]*field*q_rule.weights[iq]*abs(elt.jacob_det))
                    
                #psi=np.dot(alpha_node[i,:],w0)
                 #fvec[inode,:]=fvec[inode,:]+psi*field*q_rule.weights[iq]*abs(elt.jacob_det)
            elif tdof=="edge":
                w1=ref.hcurl_mapping(elt,
                    ref.w1_interpolation(mesh,ielt,q_rule.pts[iq,:]))
                for i in range(mesh.np_elt2n_edges[ielt]):
                    iedge=mesh.np_elt2edges[ielt,i]
                    wedge=mesh.np_elt2w_edges[ielt,i]
                    fvec[iedge]=(fvec[iedge]+
                        np.dot(w1[:,i],field)*wedge
                        *q_rule.weights[iq]*abs(elt.jacob_det))
            elif tdof=="facet":
                w2=ref.hdiv_mapping(elt,
                    ref.w2_interpolation(mesh,ielt,q_rule.pts[iq,:]))
                for i in range(mesh.np_elt2n_facets[ielt]):
                    ifacet=mesh.np_elt2facets[ielt,i]
                    wfacet=mesh.np_elt2w_facets[ielt,i]
                    fvec[ifacet]=(fvec[ifacet]+
                        np.dot(w2[:,i],field)*wfacet
                        *q_rule.weights[iq]*abs(elt.jacob_det))
            elif tdof=="element":
                w3=ref.vol_mapping(elt,
                    ref.w3_interpolation(mesh,ielt,q_rule.pts[iq,:]))
                fvec[ielt]=(fvec[ielt]+np.dot(w3[0,0],field)
                        *q_rule.weights[iq]*abs(elt.jacob_det))
        # vol control
        vol=vol+abs(elt.jacob_det)/6.

    cmat=scipy.sparse.coo_matrix(
            (cmat_val,(cmat_row,cmat_col)), shape=(n_dofs,n_dofs))

    print(io.done)    
    print(io.bcolors.INFOBLUE + '\t integrated volume: ' + str(vol) + io.bcolors.ENDC)
    
    print(' -- Solving linear systems')
    (xvec,solve_flag)=scipy.sparse.linalg.minres(cmat, fvec, tol=1e-9, maxiter=10000)
    if solve_flag==0 :
        print(io.bcolors.OKGREEN + '  >  successful ' + io.bcolors.ENDC)
    elif solve_flag>0 :
        print(io.bcolors.WARNING + '  >  max iteration # reached ' + io.bcolors.ENDC)
    elif solve_flag<0:
        print(io.bcolors.FAIL + '  >  fail ' + io.bcolors.ENDC)

    
    ## post    
    visuF=open("m2m_field.pos", 'w')
    visuF.write('General.FastRedraw = 0 ;\n')
    visuF.write('General.Color.Background = White ;\n')
    visuF.write('General.Color.Foreground = White ;\n')
    visuF.write('General.Color.Text = Black ;\n')
    visuF.write("View \"%s_field\" {\n" %('projected'))    
    #xx=np.empty(mesh.n_elements,dtype=float)
    for ielt in range(mesh.n_elements):
        elt=mesh.element_list[ielt]
        if elt.elt_type!=4: continue
        q_rule=qr.lin_tetra_1
        for iq in range(q_rule.n_pts):
            real_pos=ref.get_real_pos(elt,q_rule.pts[iq,:])
            field=np.array([0.,0.,0.])
            if tdof=="node":
                w0=ref.hgrad_mapping(elt,
                        ref.w0_interpolation(mesh,ielt,q_rule.pts[iq,:]))    
                for i in range(mesh.np_elt2n_nodes[ielt]):
                    inode=mesh.np_elt2nodes[ielt,i]
                    field[0]=field[0]+xvec[inode]*w0[0,i]
            elif tdof=="edge":
                w1=ref.hcurl_mapping(elt,
                        ref.w1_interpolation(mesh,ielt,q_rule.pts[iq,:]))    
                for i in range(mesh.np_elt2n_edges[ielt]):
                    iedge=mesh.np_elt2edges[ielt,i]
                    wedge=mesh.np_elt2w_edges[ielt,i]
                    field=field+xvec[iedge]*w1[:,i]*wedge
            elif tdof=="facet":
                w2=ref.hdiv_mapping(elt,
                        ref.w2_interpolation(mesh,ielt,q_rule.pts[iq,:]))    
                for i in range(mesh.np_elt2n_facets[ielt]):
                    ifacet=mesh.np_elt2facets[ielt,i]
                    wfacet=mesh.np_elt2w_facets[ielt,i]
                    field=field+xvec[ifacet]*w2[:,i]*wfacet
            elif tdof=="element":
                w3=ref.vol_mapping(elt,
                        ref.w3_interpolation(mesh,ielt,q_rule.pts[iq,:]))    
                field[0]=field[0]+xvec[ielt]*w3[0,0]
            visuF.write('VP(%13.5E,%13.5E,%13.5E)\n' %(real_pos[0],real_pos[1],real_pos[2]))
            visuF.write('{%13.5E,%13.5E,%13.5E};\n' %(field[0],field[1],field[2]))    
            #xx[ielt]=field[0]
    visuF.write('};\n')
    visuF.close()

    np.savetxt("m2m_field.csv",xvec)                                                
    #np.savetxt("xx.csv",xx)                                                



if __name__ == "__main__":
    main()
