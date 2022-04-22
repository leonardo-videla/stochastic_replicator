# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 22:47:14 2022

@author: leona
"""


import numpy as np
from matplotlib import pyplot as plt
import ternary
import random 

import Nabla 

def potencial (x):
    A= np.array(([0, 1, 0], [0, 0, 2], [1.5, 0, 0]))
    #B= np.dot(np.transpose (A), A)
    return np.dot (x, A.dot(x)) 
 
def potencialf1 (x):
    return np.log (1+x[0]*x[1]*x[2])
 
def potencialf2 (x):
    return 5 / (1+ 3*x[0]+2*x[1]+2*x[2])
     
def inter_potential (x, y):
    return (x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2 
    

def Replicator (x, niter):
    data=[]
    data.append(x)
    d=x
    dt=0.001
    T=dt*niter
    for i in range(niter-1):
        d=d+Nabla.CoNabla (potencial, d)*dt + 7*Nabla.CoNabla (potencialf1, d)*np.sqrt (dt)*np.random.normal (1)+ 2*Nabla.CoNabla (potencialf2, d)*np.sqrt (dt)*np.random.normal (1)
        data.append (d) 
    
    a=np.transpose(np.array(data))
    
    plt.plot (np.linspace(0, T, niter), a[0][:])
   
    
    
    plt.plot (np.linspace(0, T, niter), a[:][1])
  
    
    
    plt.plot (np.linspace(0, T, niter), a[:][2])
    plt.title("Dinámicas ")
    plt.xlabel("Tiempo")
    plt.ylabel("Abundancia")
    plt.show()
    
def Inter_Replicator (x, niter, n_replicators):
    data=[]
    dinit=[]
    data0=[]
    for i in range(n_replicators):
        x1=random.uniform (0.1, 0.9)
        x2=random.uniform (0, 0.9-x1)
        x3=1-x2-x1
       # print ("Condiciones iniciales: ")
       # print ([x1, x2, x3])
        dinit.append (np.array([x1, x2, x3]))
        if i == 0:
            data0.append (np.array([x1, x2, x3]))
            
    data.append(dinit)
        
        
    dt=0.001
    
    
    for i in range(niter-1):
        daux=dinit
        
        for j in range(n_replicators):
            dinit[j]=dinit[j]+Nabla.CoNabla (potencial, dinit[j])*dt + 7*Nabla.CoNabla (potencialf1, dinit[j])*np.sqrt (dt)*np.random.normal (1)#+ 2*Nabla.CoNabla (potencialf2, dinit[j])*np.sqrt (dt)*np.random.normal (1)
            r=np.sqrt(dt)*np.random.normal (1) 
            for k in range (n_replicators):
                dinit[j]= dinit[j] + (1/n_replicators) * 0.2*Nabla.CoNabla (inter_potential, daux[j], daux[k])*r
            if j==0:
                data0.append(dinit[j])
        
        data.append (list(dinit))
       
        #print("Datos :", dinit)
      
    
    #print (data0)
    
    a=np.transpose (np.array(data0))    
   
    return data,a
    
    #plt.plot (np.linspace(0, T, niter), a[0][:])
   
    
    
    #plt.plot (np.linspace(0, T, niter), a[:][1])
  
    
    
    #plt.plot (np.linspace(0, T, niter), a[:][2])
    #plt.title("Dinámicas ")
    #plt.xlabel("Tiempo")
    #plt.ylabel("Abundancia")
    #plt.show()
    
    
def Plot_Inter_Replicators (x, niter, n_replicators, NSAMPLES):
    
    bins=np.array([[0]*30]*30)
    for i in range (NSAMPLES):
        dat, res=Inter_Replicator (x, niter, n_replicators)
        result = [res[0][len(res[0])-1], res[1][len(res[1])-1], res[2][len(res[2])-1]]    
        j=int(np.floor(result[0]*30)-1)
        k=int(np.floor(result[1]*30)-1)
        bins[j][k]+=1
        
    def distr (p):
        j=int(np.floor(p[0]*30)-1)
        k=int(np.floor(p[1]*30)-1)
        return bins[j][k]
    
    scale=60
    figure, tax = ternary.figure(scale=scale)
    figure.set_size_inches(10, 8)
    tax.heatmapf(distr, boundary=True, style="triangular")
    tax.boundary(linewidth=2.0)
    tax.set_title("Distribucion de las abundancias relativas")
    tax.ticks(axis='lbr', linewidth=1, multiple=5)
    tax.clear_matplotlib_ticks()
    tax.show()
    
def Plot_Trajectory (x, niter, n_replicators):
        
        dat,res = Inter_Replicator (x, niter, n_replicators)
      

        figure, tax = ternary.figure(scale=1.0)
        figure.set_size_inches(5, 5)

        tax.boundary()
        tax.gridlines(multiple=0.2, color="black")
        tax.set_title("Interacting replicators", fontsize=8)
        
        
        for i in range(n_replicators):
            d=[]
            for j in range (niter):
                d.append(list(dat[j][i]))
          
            tax.plot(d, linewidth=1.5)
        
        tax.ticks(axis='lbr', multiple=0.2, linewidth=1, tick_formats="%.1f", offset=0.02)

        tax.get_axes().axis('off')
        tax.clear_matplotlib_ticks()
        tax.legend()
        tax.show()
    