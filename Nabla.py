# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 21:50:51 2022

@author: leona
"""

import numpy as np

def Nabla (fn , x, y=[]):
    n=len (x)
    h=0.001
    unit=[0]*n
    g=[]
    for i in range (n):
        unit[i]=h
        if len (y)==0:
            g.append((fn(x+np.array(unit))-fn(x))/h)
        else:
            g.append((fn(x+np.array(unit), y)-fn(x, y))/h)
        unit[i]=0
    return np.array(g)
   

def CoNabla (fn, x, y=[]):
    n=len (x)
    if len (y)==0:
        return np.multiply (x, Nabla(fn, x)-np.array([np.dot(x, Nabla(fn, x))]*n))
    else:
        return np.multiply (x, Nabla(fn, x, y)-np.array([np.dot(x, Nabla(fn, x, y))]*n))

    