#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       math_tools.py
#       
#       Copyright (C) 2013-2015 Zifu Wang <z@mesh2mesh.com>
#       


import math
import numpy as np

def cantor_pair(this):
    [x,y]=this
    return int(((x + y)*(x + y + 1)/2) + y)

def inv_cantor_pair(z):
    w=math.floor(((8*z+1)**0.5-1)/2.)
    t=(w**2+w)/2
    y=int(z-t)
    x=int(w-y)
    return [x,y]

def cantor_3(this):
    [x,y,z]=this
    a=cantor_pair([x,y])
    b=cantor_pair([a,z])
    return b

def inv_cantor_3(b):
    [a,z]=inv_cantor_pair(b)
    [x,y]=inv_cantor_pair(a)
    return [x,y,z]

def szudzik_pair(this):
    [x,y]=this
    if x<y : return y**2+x
    else : return x**2+x+y

def inv_szudzik_pair(z):
    i=int(math.floor(z**0.5))
    j=int(z-i**2)
    if i>j : return [j,i]
    else :   return [i,j-i]

def szudzik_3(this):
    [x,y,z]=this
    a=szudzik_pair([x,y])
    b=szudzik_pair([a,z])
    return b

def inv_szudzik_3(b):
    [a,z]=inv_szudzik_pair(b)
    [x,y]=inv_szudzik_pair(a)
    return [x,y,z]

def normalize_vec(x):
    y=x/(np.dot(x,x)**0.5)
    return y
