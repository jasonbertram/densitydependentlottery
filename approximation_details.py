# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:37:21 2018

@author: jbertram
"""

import numpy as np
import bisect
import matplotlib.pyplot as plt

def deltnplusnaive(m,c,U):
    if sum(m)>0 and U>0:
        L=sum(m)/float(U)
        l=m/float(U)
        cbar=sum(m*c)/sum(m)
        dr=m*np.exp(-l)*(1-np.exp(-(L-l)))*c/((cbar*L-c*l)/(1-np.exp(-(L-l)))+c)
        da=m*c*(1-np.exp(-l))/(cbar*L+c*l**2*np.exp(-l)/(1-(1+l)*np.exp(-l)))
        return m*np.exp(-L) + dr + da
    else:
        return np.zeros(len(m))

def wins(c):
    winsnocomp=np.zeros(len(c)); wins1=np.zeros(len(c)); wins2=np.zeros(len(c));
    comp=np.zeros(T);
    for i in range(int(T)):
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
                
    return winsnocomp+wins1+wins2
            


T=1000000
cs=np.array([[10.**i,1.] for i in np.linspace(-2,2,101)])

fig1, (ax1, ax2) = plt.subplots(nrows=2,figsize=[5,8])

l1=3.; l2=3.; L=l1+l2
scatter=np.random.poisson(lam=[l1,l2],size=[T,2])
m=np.array([l1*T,l2*T])

#deltaovermsim1=np.array(map(wins,cs))/m
#deltaovermnaive1=np.array(map(lambda x: deltnplusnaive(m,x,T),cs))/m
#deltaovermapprox1=np.array(map(lambda x: deltnplus(m,x,T),cs))/m

ax1.plot(deltaovermsim1[:,0],'k.',markersize=3,label=r'Simulation')
ax1.plot(deltaovermsim1[:,1],'k.',markersize=3)
ax1.plot(deltaovermnaive1[:,0],'r',linewidth=2,label=r'$\tilde{p}$,$\hat{p}$ approximation')
ax1.plot(deltaovermnaive1[:,1],'r',linewidth=2)
ax1.plot(deltaovermapprox1[:,0],'b',linewidth=2,label=r'$\tilde{q}$,$\hat{q}$ approximation')
ax1.plot(deltaovermapprox1[:,1],'b',linewidth=2)
ax1.plot([50,50],[0,1],'k',linewidth=2)
ax1.set_xticklabels([])
ax1.set_ylim([0,0.4])
ax1.set_yticklabels([0,'',0.1,'',0.2,'',0.3,'',0.4])
ax1.legend(loc='upper center',prop={'size':11})
ax1.set_ylabel(r"$\Delta_+ n_i/m_i$",fontsize=14)
ax1.annotate(r'$l_1=3$',xy=(0.1,0.12),xycoords='axes fraction',fontsize=14)
ax1.annotate(r'$l_2=3$',xy=(0.1,0.8),xycoords='axes fraction',fontsize=14)

l1=.1; l2=3.; L=l1+l2
scatter=np.random.poisson(lam=[l1,l2],size=[T,2])
m=np.array([l1*T,l2*T])

#deltaovermsim2=np.array(map(wins,cs))/m
#deltaovermnaive2=np.array(map(lambda x: deltnplusnaive(m,x,T),cs))/m
#deltaovermapprox2=np.array(map(lambda x: deltnplus(m,x,T),cs))/m

ax2.plot(deltaovermsim2[:,0],'k.',markersize=3,label=r'Simulation $l_1=0.1, l_2=3$')
ax2.plot(deltaovermsim2[:,1],'k.',markersize=3)
ax2.plot(deltaovermnaive2[:,0],'r',linewidth=2,label=r'$\tilde{p}$,$\hat{p}$ approximation')
ax2.plot(deltaovermnaive2[:,1],'r',linewidth=2)
ax2.plot(deltaovermapprox2[:,0],'b',linewidth=2,label=r'$\tilde{q}$,$\hat{q}$ approximation')
ax2.plot(deltaovermapprox2[:,1],'b',linewidth=2)
ax2.plot([50,50],[0,1],'k',linewidth=2)
ax2.set_xticks([0,25,50,75,100])
ax2.set_xticklabels([r"$10^{-2}$",r"$10^{-1}$",r"$1$",r"$10^{1}$",r"$10^{2}$"])
#ax2.set_yticklabels([r'$0$','','','','',r"$1$"])
ax2.set_ylim([0,1.])
ax2.set_ylabel(r"$\Delta_+ n_i/m_i$",fontsize=14)
ax2.set_xlabel(r"$c_1/c_2$",fontsize=16)
ax2.annotate(r'$l_1=0.1$',xy=(0.1,0.1),xycoords='axes fraction',fontsize=14)
ax2.annotate(r'$l_2=3$',xy=(0.1,0.35),xycoords='axes fraction',fontsize=14)

ax3=fig1.add_axes([0.15,0.3,0.3,0.15])

ax3.plot(deltaovermsim2,'k.',markersize=3)
ax3.plot(deltaovermnaive2,'r',linewidth=2)
ax3.plot(deltaovermapprox2,'b',linewidth=2)
ax3.plot([50,50],[0,1],'k',linewidth=2)
ax3.set_xlim([45,55])
ax3.set_ylim([0.25,0.35])
ax3.set_xticklabels(['','','','1','',''])
ax3.set_yticklabels([])

plt.savefig('/home/jbertram/repos/densitydependentlottery/approx_details.pdf',bbox="tight")


distsolo=np.array([_ for _ in scatter if _[0]==1 and _[1]==0])
distrare=np.array([_ for _ in scatter if _[0]==1 and _[1]>=1])
distcommon=np.array([_ for _ in scatter if _[0]>1])

plt.hist(cs[0,0])


#
#print np.mean(1./np.sum(distrare,1)), len(distrare)
#print np.mean(1./np.sum(distcommon,1)), len(distcommon)


#print l2/(1-np.exp(-l2))
#print l2*(L-1+np.exp(-L))/((1-(1+L)*np.exp(-L))*l2)



    
