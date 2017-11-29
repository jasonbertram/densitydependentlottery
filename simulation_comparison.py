# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:26:08 2017

@author: jbertram
"""

#Simulation comparison figure
#================================================================
m=np.array([10000,90000])
M=float(sum(m)) 
c=np.array([1.5,1.])
cbar=sum(m*c)/sum(m)
maxl=2; Us=np.array([int(m[0]/x) for x in np.linspace(0.001,maxl,1000)]);
exact=np.array([deltnplussim(m,c,U) for U in Us]);
totalclassic=np.array([deltnplusclassic(m,c,U) for U in Us]);
totalclassicboundedm=np.array([deltnplusclassicboundedm(m,c,U) for U in Us]);

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])

#First genotype
nocomp=np.array([np.exp(-M/U) for U in Us]);
comp1=np.array([R(m,c,U)[0]*c[0]/cbar for U in Us]);
comp2=np.array([A(m,c,U)[0]*c[0]/cbar for U in Us]);

ax1.plot(float(m[0])/Us,nocomp,'k',label=r"$e^{-L}$")
ax1.plot(float(m[0])/Us,comp1,'r',label=r"$R_i \frac{c_i}{\overline{c}}$")
ax1.plot(float(m[0])/Us,comp2,'b',label=r"$A_i \frac{c_i}{\overline{c}}$")
ax1.plot(float(m[0])/Us,(nocomp+comp1+comp2),'k',label=r"$e^{-L}+R_i \frac{c_i}{\overline{c}}+A_i \frac{c_i}{\overline{c}}$",linewidth=2)

ax1.plot(float(m[0])/Us,(exact[:,0,0]+exact[:,1,0]+exact[:,2,0])/m[0],'k.',markersize=1.,label="Simulation")
ax1.plot(float(m[0])/Us,exact[:,0,0]/m[0],'k.',markersize=1)
ax1.plot(float(m[0])/Us,exact[:,1,0]/m[0],'k.',markersize=1)
ax1.plot(float(m[0])/Us,exact[:,2,0]/m[0],'k.',markersize=1)

ax1.plot(float(m[0])/Us,totalclassic[:,0]/m[0],'k--',label=r"Classic lottery",linewidth=2)
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
nocomp=np.array([np.exp(-M/U) for U in Us]);
comp1=np.array([R(m,c,U)[1]*c[1]/cbar for U in Us]);
comp2=np.array([A(m,c,U)[1]*c[1]/cbar for U in Us]);
    
ax2.plot(float(m[1])/Us,nocomp,'k',label=r"$e^{-L}$")
ax2.plot(float(m[1])/Us,comp1,'r',label=r"$R_2 c_2/\overline{c}$")
ax2.plot(float(m[1])/Us,comp2,'b',label=r"$A_2 c_2/\overline{c}$")

ax2.plot(float(m[1])/Us,(exact[:,0,1]+exact[:,1,1]+exact[:,2,1])/m[1],'k.',markersize=1,label="Simulation")
ax2.plot(float(m[1])/Us,(nocomp+comp1+comp2),'k',label="Total",linewidth=2)
ax2.plot(float(m[1])/Us,totalclassic[:,1]/m[1],'k--',label=r"Classic lottery",linewidth=2)
ax2.set_xlim([0,3*maxl])
ax2.set_ylim([0,1])
ax2.xaxis.set_label_coords(0.5, -0.09)
ax2.yaxis.set_label_coords(-0.1, 0.5)
ax2.set_xlabel(r"$l_2$",fontsize=14)

ax2.plot(float(m[1])/Us,exact[:,0,1]/m[1],'k.',markersize=1,label="Simulation")
ax2.plot(float(m[1])/Us,exact[:,1,1]/m[1],'k.',markersize=1)
ax2.plot(float(m[1])/Us,exact[:,2,1]/m[1],'k.',markersize=1)

ax2.annotate('Common type',xy=(0.6,0.32),xycoords='axes fraction',fontsize=12)

plt.tight_layout()

plt.savefig('/home/jason/repos/densitydependentlottery/simulationcomparison.pdf',bbox="tight")