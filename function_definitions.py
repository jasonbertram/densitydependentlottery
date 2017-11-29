# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:25:02 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt
import bisect

def deltnplussim(m,c,U):
    scatter=np.zeros([int(U),len(m)])
    for i in range(len(m)):
        for y in xrange(int(m[i])):
            cell=int(int(U)*np.random.rand());
            scatter[cell,i]=scatter[cell,i]+1;
    
    winsnocomp=np.zeros(len(m)); wins1=np.zeros(len(m)); wins2=np.zeros(len(m));
    comp=np.zeros(U);
    for i in range(int(U)):
        comp[i]=sum(scatter[i])
        if comp[i]>0:            
            lotterycmf=np.cumsum(np.array(scatter[i])*c)
            victor=bisect.bisect(lotterycmf,np.random.rand()*lotterycmf[-1])
            
            if scatter[i][victor]==1 and comp[i]==1:
                winsnocomp[victor]=winsnocomp[victor]+1
            elif scatter[i][victor]==1:
                wins1[victor]=wins1[victor]+1
            else: 
                wins2[victor]=wins2[victor]+1
        
    return np.array([winsnocomp,wins1,wins2])

def R(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m); 
    out=cbar*np.exp(-l)*(1-np.exp(-(L-l)))\
            /(c + (L-1+np.exp(-L))/(1-(1+L)*np.exp(-L))*(cbar*L-c*l)/(L-l))
    
    for i in range(len(out)):
        if np.isnan(out)[i]: out[i]=0
            
    return out

def A(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m)
    out=cbar*(1-np.exp(-l))\
            /((1-np.exp(-l))/(1-(1+l)*np.exp(-l))*c*l\
            +(L*(1-np.exp(-L))/(1-(1+L)*np.exp(-L))-l*(1-np.exp(-l))/(1-(1+l)*np.exp(-l)))/(L-l)*(cbar*L-c*l))
    for i in range(len(out)):
        if np.isnan(out)[i]: out[i]=0
            
    return out
    
def deltnplus(m,c,U):
    if sum(m)>0 and U>0:
        L=sum(m)/float(U)
        cbar=sum(m*c)/sum(m)
        return m*(np.exp(-L)+(R(m,c,U)+A(m,c,U))*c/cbar)
    else:
        return np.zeros(len(m))
          
def deltnplusclassic(m,c,U):
    if sum(m)>0:
        return U*m*c/sum(m*c)
    else:
        return np.zeros(len(m))
        
def deltnplusclassicboundedm(m,c,U):
    if sum(m)>0:
        return np.min(np.array([m,U*m*c/sum(m*c)]),0)
    else:
        return np.zeros(len(m))     