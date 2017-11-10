# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:39:18 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt

x,y=np.linspace(0,1,1000),np.linspace(0,1,1000)
X, Y = np.meshgrid(x,y)
alphas=np.array([[1,1.5],[1.5,1]])
U = 2*(1-alphas[0,0]*X-alphas[0,1]*Y)*X
V = (1-alphas[1,1]*Y-alphas[1,0]*X)*Y
#speed = np.sqrt(U*U + V*V)

fig1, (ax1, ax2) = plt.subplots(ncols=2,figsize=[8,4])
ax1.set_aspect(1)
ax2.set_aspect(1)
plt.tight_layout()
seedpoints=np.concatenate([[[_,0.2-_] for _ in np.arange(0.02,0.2,0.02)],\
            [[float(_)/10,1.-_*0.1] for _ in range(0,11) if _ not in [2,7]]])
ax1.streamplot(x, y, U, V,density=100,start_points=seedpoints,linewidth=1.)

z1=np.linspace(0,1,10)
z2=1-z1
ax1.plot(z1,(1-alphas[0,0]*z1)/alphas[0,1],'k',linewidth=2)
ax1.plot(z1,(1-alphas[1,0]*z1)/alphas[1,1],'k',linewidth=2)
#ax1.plot([0,1/alphas[0,0]],[1/alphas[1,1],0],'y',linewidth=2)
ax1.set_xlim([0,1])
ax1.set_ylim([0,1])

ax1.annotate(r'$\frac{dn_2}{dt}=0$',xy=(0.15,0.81),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$\frac{dn_1}{dt}=0$',xy=(0.62,0.28),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$K_{21}$',xy=(0.65,-0.07),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$K_{11}$',xy=(0.95,-0.07),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$K_{12}$',xy=(-0.13,0.65),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$K_{22}$',xy=(-0.13,.97),xycoords='axes fraction',fontsize=16)
ax1.annotate(r'$(a)$',xy=(0.9,0.93),xycoords='axes fraction',fontsize=16)

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xlabel(r"$n_1$",fontsize=20)
ax1.set_ylabel(r"$n_2$",fontsize=20)

ax2.quiver(z1[:-1], z2[:-1], z1[1:]-z1[:-1], z2[1:]-z2[:-1], scale_units='xy', angles='xy', scale=1)
ax2.annotate(r'$\frac{dN}{dt}=0$',xy=(0.53,0.53),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$(b)$',xy=(0.9,0.93),xycoords='axes fraction',fontsize=16)

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xlabel(r"$n_1$",fontsize=20)
ax2.set_ylabel(r"$n_2$",fontsize=20)

plt.savefig('/home/jbertram/repos/densitydependentlottery/Kplot.pdf',bbox="tight")



#fig1, ax1 = plt.subplots()
#ax1.streamplot(X, Y, U, V,density=0.5)
#plt.xlim([0,2])
#plt.ylim([0,2])
#
#x=np.linspace(0,2,100)
#plt.plot(x,1-0.5*x,'k',linewidth=2)
#plt.plot(x,2-2*x,'k',linewidth=2)