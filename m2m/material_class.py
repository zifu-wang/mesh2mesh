#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       material_class.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#       
import os
import sys
import numpy as np
import m2m.io_base as io

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.integrate import quad

# The absolute path of the directoy for this module file:
_ROOT = os.path.abspath(os.path.dirname(__file__))

mu_0=4*np.pi*1e-7

class Material(object):
    def __init__(self, name=None, kind=None, tag=None,
            mu_r=None, bhfile=None, 
            sigma=None):
        self.name = name
        self.kind = kind    # linear/nonlinear
        self.tag = tag      # corresponded num in element tags
        if kind =='linear': # linear mat has fixed permeability 
            self.mu_r = mu_r
        elif kind =='nonlinear':  
            bhtab = np.loadtxt(bhfile) # B-H curve file needed
            self.hdata = bhtab[:,0]
            self.bdata = bhtab[:,1]
            # B-H interpolation
            self.interp_h2b = interp1d(self.hdata,self.bdata,'quadratic')
            self.interp_b2h = interp1d(self.bdata,self.hdata,'quadratic')
            # B-H linear regression
            lpopt,lpopc = curve_fit(linear_h2b,self.hdata,self.bdata)
            self.mu_r = lpopt[0]
            # B-H marrocco curve fit
            nlpopt,nlpopc = curve_fit(marrocco_b2h,self.bdata,self.hdata,
                    p0=[1e+1, 5e+6, 2e-1, 2e-3])
            self.marrocco_values = nlpopt
        else: sys.exit(io.bcolors.FAIL + 'Wrong material kind ' + io.bcolors.ENDC) 
        self.sigma = sigma
    def h2b(self,h):
        if self.kind == 'linear':
            return linear_h2b(h, self.mu_r)
        elif self.kind == 'nonlinear':
            h_norm = np.dot(h,h)**0.5
            try:    return h/h_norm*(self.interp_h2b(h_norm))
            except:
                (alpha,tau,c,eps) = self.marrocco_values
                b_norm = marrocco_h2b(h_norm,alpha,tau,c,eps)
                return h/h_norm*b_norm
    def b2h(self,b):
        if self.kind == 'linear':
            return linear_b2h(b, self.mu_r)
        elif self.kind == 'nonlinear':
            b_norm = np.dot(b,b)**0.5
            try:     return b/b_norm*(self.interp_b2h(b_norm))
            except: 
                (alpha,tau,c,eps) = self.marrocco_values
                h_norm = marrocco_b2h(b_norm,alpha,tau,c,eps)
                return b/b_norm*h_norm
    def b2coenergy(self,b):
        b_norm=np.dot(b,b)**0.5
        if self.kind == 'linear':
            return 0.5*b_norm/(mu_0*self.mu_r)
        elif self.kind == 'nonlinear':
            (coenergy, error) = quad (self.b2h, 0, b_norm, epsrel=1e-5)
            return coenergy
    def b2energy(self,b):
        b_norm=np.dot(b,b)**0.5
        if self.kind == 'linear':
            return 0.5*b_norm/(mu_0*self.mu_r)
        elif self.kind == 'nonlinear':
            (energy, error) = quad (self.h2b, 0, self.b2h(b_norm), epsrel=1e-5)
            return energy
        
def linear_h2b(h,mu_r):
    return mu_r*mu_0*h

def linear_b2h(b,mu_r):
    return b/mu_r/mu_0

def marrocco_b2h(b,alpha,tau,c,eps):
    return b/mu_0*(b**(2*alpha)/(b**(2*alpha)+tau)*(c-eps)+eps)

def marrocco_h2b(h,alpha,tau,c,eps,b_max=10e0,b_min=0e0,tol=1e-3,max_ite=99):
    uboundary=b_max
    lboundary=b_min
    for i in range(max_ite):
        b_test=(uboundary+lboundary)/2
        h_test=marrocco_b2h(b_test,alpha,tau,c,eps)
        if abs(h_test-h)/abs(h)<=tol: return b_test
        elif h_test<h: lboundary=b_test
        elif h_test>h: uboundary=b_test
    return (uboundary+lboundary)/2

# Material library
air=Material(    
        name='air',
        kind='linear',
        mu_r=1.0,
        )

core=Material(    
        name='core',
        kind='linear',
        mu_r=1e3,
        )

copper=Material(    
        name='copper',
        kind='linear',
        mu_r=0.99999e0,
        sigma=5.8e7,
        )

m15=Material(    
        name='m15',
        kind='nonlinear',
        bhfile=os.path.join(_ROOT, 'materials/m15_bh.tab'),
        )

m27=Material(    
        name='m27',
        kind='nonlinear',
        bhfile=os.path.join(_ROOT, 'materials/m27_bh.tab'),
        )

steel1010=Material(    
        name='steel1010',
        kind='nonlinear',
        bhfile=os.path.join(_ROOT, 'materials/steel1010_bh.tab'),
        )


