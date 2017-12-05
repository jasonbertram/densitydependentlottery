# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:29:53 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#SS figure
#========================================================

fig1,ax1=plt.subplots(figsize=[4,4])

s=np.linspace(0,0.5,100)
ax1.plot(s,1/(1-s),'k',linewidth=2)

ax1.set_ylim([1,2])
ax1.set_xlabel(r"$\epsilon$",fontsize=20)
ax1.set_ylabel(r"$s_{\rm final}/s_{\rm initial}$",fontsize=20)

plt.tight_layout()


plt.savefig('/home/jbertram/repos/densitydependentlottery/strengthofselection.pdf',bbox="tight")

#Logistic model
#==========================================================

r=[1,1];
s=0.3
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
totaltime=80
b=np.array([3.,3.])
c=np.array([1.,1.])
d=np.array([1.,1.])
n0=np.array([49000.,100.])

n=n0;
nhist=[n];
for t in range(totaltime/2):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))

#plt.plot(np.sum(nhist,1))
n0=np.array(n)
print n0 #equilibrium

s=0.4
d=np.array([1.,1.*(1-s)])
n=n0;
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))
    
nhist=np.array(nhist)
Nhist=np.sum(nhist,1)
print "b absolute: ", nhist[1,1]/nhist[0,1]-1
print "b relative: ", nhist[1,1]/nhist[0,1]*Nhist[0]/Nhist[1]-1
print "N change: ", Nhist[-1]/Nhist[0]


fig1,ax1=plt.subplots()
ax1.plot(nhist[:,1]/Nhist)

fig2,ax2=plt.subplots()
ax2.plot(Nhist)


d=np.array([1.,1.])
c=np.array([1.,3.])
n=n0;
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))
    
nhist=np.array(nhist)
Nhist=np.sum(nhist,1)
print nhist[1,1]/nhist[0,1]-1
print nhist[1,1]/nhist[0,1]*Nhist[0]/Nhist[1]-1

ax1.plot(nhist[:,1]/Nhist)
ax2.plot(Nhist)



#plt.figure()
#plt.plot(map(np.log,nhist[:,1]))

#plt.figure()
#plt.plot(np.sum(nhist,1))
#plt.ylim([0,T])