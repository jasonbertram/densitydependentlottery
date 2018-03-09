# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:26:08 2017

@author: jbertram
"""

#Simulation comparison figure
#================================================================
c=np.array([1.5,1.])
cbar=sum(np.array([0.1,0.9])*c)
U=100000
ls=np.array([[x*0.1,x] for x in np.linspace(0.01,20,1000)])

exact=np.array([deltnplussim(m,c,U) for m in U*ls]);
totalclassic=np.array([deltnplusclassic(m,c,U) for m in U*ls]);
#totalclassicboundedm=np.array([deltnplusclassicboundedm(m,c,U) for m in U*ls]);

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

#First genotype
nocomp=np.array([np.exp(-np.sum(m)/U) for m in U*ls]);
comp1=np.array([R(m,c,U)[0]*c[0]/cbar for m in U*ls]);
comp2=np.array([A(m,c,U)[0]*c[0]/cbar for m in U*ls]);

ax1.plot(ls[:,0],nocomp,'k',label=r"$e^{-L}$")
ax1.plot(ls[:,0],comp1,'r',label=r"$R_i \frac{c_i}{\overline{c}}$")
ax1.plot(ls[:,0],comp2,'b',label=r"$A_i \frac{c_i}{\overline{c}}$")
ax1.plot(ls[:,0],(nocomp+comp1+comp2),'k',label=r"$e^{-L}+R_i \frac{c_i}{\overline{c}}+A_i \frac{c_i}{\overline{c}}$",linewidth=2)

ax1.plot(ls[:,0],(exact[:,0,0]+exact[:,1,0]+exact[:,2,0])/(ls[:,0]*U),'k.',markersize=1.,label="Simulation")
ax1.plot(ls[:,0],exact[:,0,0]/(ls[:,0]*U),'k.',markersize=1)
ax1.plot(ls[:,0],exact[:,1,0]/(ls[:,0]*U),'k.',markersize=1)
ax1.plot(ls[:,0],exact[:,2,0]/(ls[:,0]*U),'k.',markersize=1)

ax1.plot(ls[:,0],totalclassic[:,0]/(ls[:,0]*U),'k--',label=r"Classic lottery",linewidth=2)
ax1.set_xlim([0,maxl])
ax1.set_ylim([0,1])
ax1.xaxis.set_label_coords(0.5, -0.09)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xticklabels(['0','','1','','2','','3'])
ax1.set_xlabel(r"$l_1$",fontsize=14)
ax1.set_ylabel("Success per propagule",fontsize=12)

ax1.annotate('Rare type',xy=(0.65,0.32),xycoords='axes fraction',fontsize=12)
ax1.legend(loc='upper center',prop={'size':11})

#Second genotype
nocomp=np.array([np.exp(-np.sum(m)/U) for m in U*ls]);
comp1=np.array([R(m,c,U)[1]*c[1]/cbar for m in U*ls]);
comp2=np.array([A(m,c,U)[1]*c[1]/cbar for m in U*ls]);
    
ax2.plot(ls[:,1],nocomp,'k',label=r"$e^{-L}$")
ax2.plot(ls[:,1],comp1,'r',label=r"$R_2 c_2/\overline{c}$")
ax2.plot(ls[:,1],comp2,'b',label=r"$A_2 c_2/\overline{c}$")

ax2.plot(ls[:,1],(exact[:,0,1]+exact[:,1,1]+exact[:,2,1])/(ls[:,1]*U),'k.',markersize=1,label="Simulation")
ax2.plot(ls[:,1],(nocomp+comp1+comp2),'k',label="Total",linewidth=2)
ax2.plot(ls[:,1],totalclassic[:,1]/(ls[:,1]*U),'k--',label=r"Classic lottery",linewidth=2)
ax2.set_xlim([0,20])
ax2.set_ylim([0,1])
ax2.xaxis.set_label_coords(0.5, -0.09)
ax2.yaxis.set_label_coords(-0.1, 0.5)
ax2.set_xlabel(r"$l_2$",fontsize=14)

ax2.plot(ls[:,1],exact[:,0,1]/(ls[:,1]*U),'k.',markersize=1,label="Simulation")
ax2.plot(ls[:,1],exact[:,1,1]/(ls[:,1]*U),'k.',markersize=1)
ax2.plot(ls[:,1],exact[:,2,1]/(ls[:,1]*U),'k.',markersize=1)

ax2.annotate('Common type',xy=(0.6,0.32),xycoords='axes fraction',fontsize=12)

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/simulationcomparison.pdf',bbox="tight")