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

fig1, ax1 = plt.subplots(ncols=1,figsize=[4,4])

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

#def f(n,t):
#    return (b-d*sum(n))*n
#    
#def g(p,t):
#    return s*p*(1-p)
#
#solK=odeint(f,[0.99*b[0]/d[0],0.01*b[0]/d[0]],x)
#solrel=odeint(g,0.01,x)
#
#ax2.plot(x,solK[:,1]/(solK[:,0]+solK[:,1]),'k--', linewidth=2,label=r"$\epsilon=0.2$ (Eq. 11)")
#ax2.plot(x,solrel,'k', linewidth=2,label=r"Canonical (Eq. 1)")
#ax2.legend(loc='lower right',prop={'size':11})
#ax2.annotate(r'$(b)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
#ax2.set_xticklabels([])
#ax2.set_xlabel("Time",fontsize=14)
#ax2.set_yticklabels(['0','','','','','1'])
#ax2.set_ylabel("Frequency",fontsize=14)

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/strengthofselection.pdf',bbox="tight")


#lottery bsweep
#==========================================================

T=100000
totaltime=100
b=np.array([1.,1.])
c=np.array([1.,1.])
d=np.array([1.5,1.5])
n0=np.array([1000.,1000.])

n=n0;
nhist=[n];
for t in range(totaltime/2):
    U=T-sum(n)
    n=n*(1+deltnplus(b*n*U/T,c,U)/n)/d
    nhist.append(list(n))

print nhist[-1]

totaltime=100
b=np.array([1.,1.5])
c=np.array([1.,1.])
d=np.array([1.5,1.5])
n0=np.array([sum(nhist[-1]),10.])

n=n0;
nhist=[n]; Whist=[]
for t in range(totaltime):
    U=T-sum(n)
    Whist.append((1+deltnplus(b*n*U/T,c,U)/n)/d)
    n=n*Whist[-1]
    nhist.append(list(n))

nhist=np.array(nhist)
Nhist=np.sum(nhist,1)
Whist=np.array(Whist)

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

ax1.plot(nhist[:,1],"k", linewidth=2,label=r"$b_2=1.5$")
ax1.plot(nhist[:,0],"k--", linewidth=2,label=r"$b_1=1$")
ax1.plot(Nhist,"k:", linewidth=2,label=r"Total")
ax1.annotate(r'$(a)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax1.set_xticklabels([])
ax1.set_xlabel("Time",fontsize=14)
ax1.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
ax1.set_ylabel(r"Density/$T$",fontsize=14)
ax1.set_ylim([0,70000])
ax1.legend(loc='upper center',prop={'size':11})

Wbarhist=np.sum((nhist[:-1]*Whist),1)/Nhist[:-1]
ax2.plot((Whist[:,1]-Whist[:,0])/Wbarhist,'k',linewidth=2,label=r"Density-dependent")
s=b[1]/b[0]-1
f=d[0]-1
ax2.plot(f/(1+f)*s/(1+s*nhist[:-1,1]/Nhist[:-1]),'k--',linewidth=2,label=r"$f(\overline{b},N)=\mathrm{constant}$")
ax2.annotate(r'$(b)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax2.set_ylim([0.1,0.19])
ax2.set_xticklabels([])
ax2.set_xlabel("Time",fontsize=14)
#ax2.set_ylim([0.2,0.25])
#ax2.set_yticklabels(['0','','','','','1'])
ax2.set_ylabel(r"$(W_2-W_1)/\overline{W}$")
ax2.legend(loc='upper center',prop={'size':11})

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/bsweep.pdf',bbox="tight")

#lottery multiple traits
#==========================================================

T=100000
totaltime=100
b=np.array([1.,1.])
c=np.array([1.,1.])
d=np.array([1.5,1.5])
n0=np.array([1000.,1000.])

n=n0;
nhist=[n];
for t in range(totaltime/2):
    U=T-sum(n)
    n=n*(1+deltnplus(b*n*U/T,c,U)/n)/d
    nhist.append(list(n))

print nhist[-1]

totaltime=100
b=np.array([1.,1.25])
c=np.array([1.,1.])
d=np.array([1.5,1.4])
n0=np.array([sum(nhist[-1]),10.])

n=n0;
nhist=[n]; Whist=[]
for t in range(totaltime):
    U=T-sum(n)
    Whist.append((1+deltnplus(b*n*U/T,c,U)/n)/d)
    n=n*Whist[-1]
    nhist.append(list(n))

nhist=np.array(nhist)
Nhist=np.sum(nhist,1)
Whist=np.array(Whist)

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

ax1.plot(nhist[:,1],"k", linewidth=2,label=r"$b_2=1.25$, $d_2=1.4$")
ax1.plot(nhist[:,0],"k--", linewidth=2,label=r"$b_1=1$, $d_1=1.5$")
ax1.plot(Nhist,"k:", linewidth=2,label=r"Total")
ax1.annotate(r'$(a)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax1.set_xticklabels([])
ax1.set_xlabel("Time",fontsize=14)
ax1.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
ax1.set_ylabel(r"Density/$T$",fontsize=14)
ax1.set_ylim([0,80000])
ax1.legend(loc='upper center',prop={'size':11})

Wbarhist=np.sum((nhist[:-1]*Whist),1)/Nhist[:-1]
ax2.plot((Whist[:,1]-Whist[:,0])/Wbarhist,'k',linewidth=2,label=r"Density-dependent")

bbarhist=np.sum((nhist[:-1]*b),1)/Nhist[:-1]
f=(1-np.exp(-bbarhist[0]*Nhist[0]/T))*(T-Nhist[0])/Nhist[0]
W1histfconst=(1+f*b[0]/bbarhist)/d[0]
W2histfconst=(1+f*b[1]/bbarhist)/d[1]
Wbarhistfconst=(W1histfconst*nhist[:-1,0]+W1histfconst*nhist[:-1,1])/Nhist[:-1]
ax2.plot((W2histfconst-W1histfconst)/Wbarhistfconst,'k--',linewidth=2,label=r"$f(\overline{b},N)=\mathrm{constant}$")
ax2.annotate(r'$(b)$',xy=(0.01,0.93),xycoords='axes fraction',fontsize=16)
ax2.set_ylim([0.11,0.18])
ax2.set_xticklabels([])
ax2.set_xlabel("Time",fontsize=14)
#ax2.set_ylim([0.2,0.25])
#ax2.set_yticklabels(['0','','','','','1'])
ax2.set_ylabel(r"$(W_2-W_1)/\overline{W}$")
ax2.legend(loc='upper center',prop={'size':11})

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/bdsweep.pdf',bbox="tight")


