# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:40:51 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

T=100000
totaltime=80
b=np.array([1.,1.])
c=np.array([1.,1.])
d=np.array([0.1,0.1])
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

s=1.
b=np.array([1.,1.*(1+s)])
n=n0;
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))
    
nhist=np.array(nhist)
Nhist=np.sum(nhist,1)
Dhist=Nhist/T
bhist=np.array(np.sum(b*nhist,1)/Nhist)
bhist[0]=1.

plt.plot((1/Dhist-1)*(1-np.exp(-bhist*Dhist))/d[0])

plt.plot(Dhist)

#=======================================================
fig1, ax1 = plt.subplots(figsize=[5,4])

T=100000
b=np.array([1.,1.])
c=np.array([1.,1.])
d=np.array([.5,.5])
eps=0.1

ns=np.array([_*np.array([1,1]) for _ in np.linspace(1,T/2,100)]) 

bWs=np.array([deltnplus(np.array([1.,1+eps])*b*x*(1-sum(x)/T),c,T-sum(x))/x-d for x in ns])
bDDS=np.sum(np.array([-1,1])*bWs,1)/(1+np.mean(bWs,1))
cWs=np.array([deltnplus(b*x*(1-sum(x)/T),np.array([1.,1+eps])*c,T-sum(x))/x-d for x in ns])
cDDS=np.sum(np.array([-1,1])*cWs,1)/(1+np.mean(cWs,1))
dWs=np.array([deltnplus(b*x*(1-sum(x)/T),c,T-sum(x))/x-np.array([1.+eps,1.])*d for x in ns])
dDDS=np.sum(np.array([-1,1])*dWs,1)/(1+np.mean(dWs,1))
ax1.plot(np.sum(ns,1)/T,bDDS,'k',linewidth=2,label=r"$b$ selection")
ax1.plot(np.sum(ns,1)/T,cDDS,'k--',linewidth=2,label=r"$c$ selection")
ax1.plot(np.sum(ns,1)/T,dDDS,'k:',linewidth=2,label=r"$d$ selection")
ax1.legend(loc='upper center',prop={'size':11})
#ax1.set_ylim([0,0.1])

#eq_dens=optimize.fsolve(lambda D: (1-np.exp(-b[0]*D))*(1/D-1)-d[0],0.5)[0]
#ax1.plot([eq_dens,eq_dens],[0,0.1],'k',linewidth=2)

ax1.set_xlabel(r"Density $N/T$")
ax1.set_ylabel(r"$\Delta n_2/n_2-\Delta n_1/n_1$")

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/DDS_lottery.pdf',bbox="tight")


