#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       field_projection.py
#       
#       Copyright 2013 Zifu Wang <zifu.wang@icloud.com>
#

import sys
#sys.setrecursionlimit(10000) # enable if kdtree runtime_error

import numpy as np
import io_base as io

import mesh_tree
import quadrature_class as qr
import reference_mapping as ref
import element_matrix as eltm

import scipy.sparse
import scipy.sparse.linalg

def main(smesh,sfield,mesh,tdof=None):
	""" Mesh-to-Mesh field projection """
	""" Note: meshes should be pre_processed before calling """
	print io.line
	print (io.bcolors.HEADER + '    Projecting field on: ' + tdof + ' dofs\n' 
							      + io.bcolors.ENDC)

        print ' -- Preparing memory'
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
	print io.done

	print ' -- Locating sample points'
	# first loop, count number
	n_smpl_pts=0
	for ielt in xrange(mesh.n_elements):
		elt=mesh.element_list[ielt]
		if elt.elt_type!=4: continue
		q_rule=qr.lin_tetra_2
		n_smpl_pts=n_smpl_pts+q_rule.n_pts
	pts=np.empty([n_smpl_pts,3],dtype=float)
	# second loop, stock pts
	i_pt=0
	for ielt in xrange(mesh.n_elements):
		elt=mesh.element_list[ielt]
		if elt.elt_type!=4: continue
		for iq in xrange(q_rule.n_pts):
			pts[i_pt,:]=(ref.get_real_pos(elt,q_rule.pts[iq,:]))
			i_pt=i_pt+1
	# now call the searching routine
	pts2elts=mesh_tree.fort_find_containing_elt(smesh,pts)

        print ' -- Assemblying matrix'
	vol=0.
	i_pt=0
	cmat_index=0
	for ielt in xrange(mesh.n_elements):
		elt=mesh.element_list[ielt]
		if elt.elt_type!=4: continue

		## Cmat
		if tdof=="node":
			cmat_grow=4**2
			(row,col,val)=eltm.w0_w0(elt)
		elif tdof=="edge":
			cmat_grow=6**2
			(row,col,val)=eltm.w1_w1(elt)
		elif tdof=="facet":
			cmat_grow=4**2
			(row,col,val)=eltm.w2_w2(elt)
		elif tdof=="element":
			cmat_grow=1
			(row,col,val)=eltm.w3_w3(elt,ielt)

		#assembly
		cmat_row[cmat_index:cmat_index+cmat_grow]=row
		cmat_col[cmat_index:cmat_index+cmat_grow]=col
		cmat_val[cmat_index:cmat_index+cmat_grow]=val
		cmat_index=cmat_index+cmat_grow

		#node biort
#		for i in xrange(elt.n_nodes):
#			inode=elt.nodes[i]
#			cmat[inode,0]=cmat[inode,0]+abs(elt.jacob_det)
		#elt
#		cmat[ielt,0]=6./abs(elt.jacob_det)
		#q_rule=quadrature_class.lin_tetra_2
		#for iq=xrange(q_rule.n_pts)		
		#	w_elt=6/abs(elt.jacob_det)
		#	cmat[ielt]=cmat[ielt]+
		#		w_elt**2*q_rule.weights[iq]*abs(elt.jacob_det)
		
		## Fvec
		q_rule=qr.lin_tetra_2
		for iq in xrange(q_rule.n_pts):
			# Quadrature pt
			#real_pos=ref.get_real_pos(elt,q_rule.pts[iq,:])
			# locate on source mesh
			iselt= pts2elts[i_pt]
			i_pt=i_pt+1
			if iselt==-1: continue
			# field
#			field=np.array([real_pos[0]**2,0,0],dtype=float)
#			field=real_pos[0]**2
			#field=[0,sfield[iselt],0]	# todo: vector field & interpolation
			field=sfield[iselt]	# todo: vector field & interpolation

			if tdof=="node":
				w0=ref.hgrad_mapping(elt,
					ref.w0_interpolation(elt,q_rule.pts[iq,:]))
				for i in xrange(elt.n_nodes):
					inode=elt.nodes[i]
				 	fvec[inode]=(fvec[inode]+
						w0[0,i]*field
						*q_rule.weights[iq]*abs(elt.jacob_det))
					
				#psi=np.dot(alpha_node[i,:],w0)
			 	#fvec[inode,:]=fvec[inode,:]+psi*field*q_rule.weights[iq]*abs(elt.jacob_det)
			elif tdof=="edge":
				w1=ref.hcurl_mapping(elt,
					ref.w1_interpolation(elt,q_rule.pts[iq,:]))
				for i in xrange(elt.n_edges):
					iedge=elt.edges[i]
					fvec[iedge]=(fvec[iedge]+
						np.dot(w1[:,i],field)*elt.w_edges[i]
						*q_rule.weights[iq]*abs(elt.jacob_det))
			elif tdof=="facet":
				w2=ref.hdiv_mapping(elt,
					ref.w2_interpolation(elt,q_rule.pts[iq,:]))
				for i in xrange(elt.n_facets):
					ifacet=elt.facets[i]
					fvec[ifacet]=(fvec[ifacet]+
						np.dot(w2[:,i],field)*elt.w_facets[i]
						*q_rule.weights[iq]*abs(elt.jacob_det))
			elif tdof=="element":
				w3=ref.vol_mapping(elt,
					ref.w3_interpolation(elt,q_rule.pts[iq,:]))
				fvec[ielt]=(fvec[ielt]+np.dot(w3[0,0],field)
						*q_rule.weights[iq]*abs(elt.jacob_det))
		# vol control
		vol=vol+abs(elt.jacob_det)/6.

	cmat=scipy.sparse.coo_matrix(
			(cmat_val,(cmat_row,cmat_col)), shape=(n_dofs,n_dofs))

	print io.done,	
	print io.bcolors.INFOBLUE + '\t integrated volume: ' + str(vol) + io.bcolors.ENDC
	
        print ' -- Solving linear systems'
	(xvec,solve_flag)=scipy.sparse.linalg.minres(cmat, fvec, tol=1e-9, maxiter=10000)
	if solve_flag==0 :
		print io.bcolors.OKGREEN + '  >  successful ' + io.bcolors.ENDC
	elif solve_flag>0 :
		print io.bcolors.WARNING + '  >  max iteration # reached ' + io.bcolors.ENDC
	elif solve_flag<0:
		print io.bcolors.FAIL + '  >  fail ' + io.bcolors.ENDC

	
	## post	
	visuF=open("test.pos", 'w')
	visuF.write('General.FastRedraw = 0 ;\n')
	visuF.write('General.Color.Background = White ;\n')
	visuF.write('General.Color.Foreground = White ;\n')
	visuF.write('General.Color.Text = Black ;\n')
	visuF.write("View \"%s_field\" {\n" %('projected'))	
	for ielt in xrange(len(mesh.element_list)):
		elt=mesh.element_list[ielt]
		if elt.elt_type!=4: continue
		q_rule=qr.lin_tetra_1
		for  iq in xrange(q_rule.n_pts):
			real_pos=ref.get_real_pos(elt,q_rule.pts[iq,:])
			field=np.array([0.,0.,0.])
			if tdof=="node":
				w0=ref.hgrad_mapping(elt,
			       		ref.w0_interpolation(elt,q_rule.pts[iq,:]))	
				for i in xrange(elt.n_nodes):
        		        	inode=elt.nodes[i]
					field[0]=field[0]+xvec[inode]*w0[0,i]
			elif tdof=="edge":
				w1=ref.hcurl_mapping(elt,
					ref.w1_interpolation(elt,q_rule.pts[iq,:]))	
				for i in xrange(elt.n_edges):
        	        	        iedge=elt.edges[i]
					field=field+xvec[iedge]*w1[:,i]*elt.w_edges[i]
			elif tdof=="facet":
				w2=ref.hdiv_mapping(elt,
					ref.w2_interpolation(elt,q_rule.pts[iq,:]))	
				for i in xrange(elt.n_facets):
        	        	        ifacet=elt.facets[i]
					field=field+xvec[ifacet]*w2[:,i]*elt.w_facets[i]
			elif tdof=="element":
				w3=ref.vol_mapping(elt,
					ref.w3_interpolation(elt,q_rule.pts[iq,:]))	
				field[0]=field[0]+xvec[ielt]*w3[0,0]
		#elt
#		w=6./abs(elt.jacob_det)
#		field=xvec[ielt]*w
			visuF.write('VP(%13.5E,%13.5E,%13.5E)\n' %(real_pos[0],real_pos[1],real_pos[2]))
			visuF.write('{%13.5E,%13.5E,%13.5E};\n' %(field[0],field[1],field[2]))	
	visuF.write('};\n')
	visuF.close()

	np.savetxt("m2m_field.csv",xvec)												



if __name__ == "__main__":
	main()
