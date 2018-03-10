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

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

s=np.linspace(0,0.5,100)
ax1.plot(s,1/(1-s),'k',linewidth=2)

ax1.set_ylim([1,2])
ax1.set_xlabel(r"$\epsilon$",fontsize=20)
ax1.set_ylabel(r"$s_{\rm final}/s_{\rm initial}$",fontsize=20)

ax1.xaxis.set_label_coords(0.5, -0.09)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.annotate(r'$(a)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)

b=np.array([1,1]);
s=0.2
d=np.array([1e-4,1e-4*(1-s)]);

x=np.linspace(0,10/s,1000)

def f(n,t):
    return (b-d*sum(n))*n
    
def g(p,t):
    return s*p*(1-p)

solK=odeint(f,[0.99*b[0]/d[0],0.01*b[0]/d[0]],x)
solrel=odeint(g,0.01,x)

ax2.plot(x,solK[:,1]/(solK[:,0]+solK[:,1]),'k--', linewidth=2,label=r"$\epsilon=0.2$ (Eq. 11)")
ax2.plot(x,solrel,'k', linewidth=2,label=r"Canonical (Eq. 1)")
ax2.legend(loc='lower right',prop={'size':11})
ax2.annotate(r'$(b)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax2.set_xticklabels([])
ax2.set_xlabel("Time",fontsize=14)
ax2.set_yticklabels(['0','','','','','1'])
ax2.set_ylabel("Frequency",fontsize=14)

plt.tight_layout()

plt.savefig('/home/jason/repos/densitydependentlottery/strengthofselection.pdf',bbox="tight")


#lottery multiple traits
#==========================================================

T=100000
totaltime=100
b=np.array([1.,1.])
c=np.array([1.,1.])
d=np.array([.9,.9])
n0=np.array([1000.,1000.])

n=n0;
nhist=[n];
for t in range(totaltime/2):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))

print nhist[-1]

totaltime=80
b=np.array([1.,1.])
c=np.array([1.,2.])
d=np.array([.9,.7])
n0=np.array([sum(nhist[-1]),10.])

n=n0;
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n=n+deltnplus(b*n*U/T,c,U)-d*n
    nhist.append(list(n))

nhist=np.array(nhist)
Nhist=np.sum(nhist,1)

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

ax1.plot(nhist[:,1],"k", linewidth=2,label=r"$d_j=0.7$, $c_j=2$")
ax1.plot(nhist[:,0],"k--", linewidth=2,label=r"$d_i=0.9$, $c_i=1$")
ax1.plot(Nhist,"k:", linewidth=2,label=r"Total")
ax1.annotate(r'$(a)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax1.set_xticklabels([])
ax1.set_xlabel("Time",fontsize=14)
ax1.set_yticklabels([0,0.05,0.1,0.15,0.2,0.25,0.3])
ax1.set_ylabel(r"Density/$T$",fontsize=14)
ax1.set_ylim([0,30000])
ax1.legend(loc='upper center',prop={'size':11})

ax2.plot((nhist[1:,1]-nhist[0:-1,1])/nhist[0:-1,1]-(nhist[1:,0]-nhist[0:-1,0])/nhist[0:-1,0],'k',linewidth=2)
#plt.ylim([0,0.5])
ax2.annotate(r'$(b)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax2.set_xticklabels([])
ax2.set_xlabel("Time",fontsize=14)
ax2.set_ylim([0.2,0.25])
#ax2.set_yticklabels(['0','','','','','1'])
ax2.set_ylabel(r"$\Delta n_j/n_j-\Delta n_i/n_i$")

plt.tight_layout()

plt.savefig('/home/jason/repos/densitydependentlottery/multiple.pdf',bbox="tight")

