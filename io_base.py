#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       io_base.py
#       
#       Copyright 2013 Zifu Wang <z@mesh2mesh.com>
#       

import os


# colorful output is always better 
class bcolors:
        HEADER = "\033[1m"
        WARNING = '\033[95m'
        INFOBLUE = '\033[94m'
        INFOYELLOW = '\033[93m'
        OKGREEN = '\033[92m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'

# try to convert a string into integer or float 
def convert(n):
        try:   return int(n)
        except ValueError:
                try:    return float(n)
                except  ValueError:	return n


line='\n-------------------------------------------------------'
done=bcolors.OKGREEN + '  >  ok' + bcolors.ENDC 
