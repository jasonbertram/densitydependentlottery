# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:29:53 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Logistic model
#==========================================================

r=[1,1];
s=0.1
K=[1e5,1e5*(1+s)];

x=np.linspace(0,20/s,1000)

def f(n,t):
    return [r[0]*(1-(n[0]+n[1])/K[0])*n[0],r[0]*(1-(n[0]+n[1])/K[1])*n[1]]
    
def g(p,t):
    return s*p*(1-p)

solK=odeint(f,[K[0]-1.,1.],x)
solrel=odeint(g,1/K[0],x)

plt.plot(x,solK[:,1]/(solK[:,0]+solK[:,1]))
plt.plot(x,solrel)
plt.ylim([0,1])

#lottery model
#==========================================================

T=100000
totaltime=200
s=0.5
b=np.array([1.,1.*(1+s)])
c=np.array([1.,1.])
d=np.array([.1,.1])

n=np.array([90000.,100.])
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))
    
nhist=np.array(nhist)

plt.plot(nhist[:,1])

plt.figure()
plt.plot(np.sum(nhist,1))

plt.ylim([0,T])

plt.figure()
plt.plot(map(np.log,nhist[:,1]))

print (nhist[11,1]-nhist[10,1])/nhist[10,1]




n=np.array([90000.,100.])
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+np.sum(deltnplussim(b*n*U/T,c,U),0)-d*n
    nhist.append(list(n))
    
nhist=np.array(nhist)

plt.plot(nhist[:,1])

plt.figure()
plt.plot(np.sum(nhist,1))