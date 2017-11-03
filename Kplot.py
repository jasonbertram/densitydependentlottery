# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:39:18 2017

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt

Y, X = np.mgrid[0:1:1000j, 0:1:1000j]
U = 2*(1-X-1.5*Y)*X
V = (1-Y-1.5*X)*Y
speed = np.sqrt(U*U + V*V)

fig1, ax1 = plt.subplots()
ax1.set_aspect(1)
#ax1.quiver(X, Y, U/speed, V/speed,alpha=0.5)
seedpoints=np.array([[1,1]])
ax1.streamplot(X, Y, U, V,color='k',start_points=)
plt.xlim([0,1])
plt.ylim([0,1])

x=np.linspace(0,2,100)
plt.plot(x,(1-x)/1.5,'k',linewidth=2)
plt.plot(x,1-1.5*x,'k',linewidth=2)

ax1.annotate(r'$f_2(n_1,n_2)=0$',xy=(0.15,0.815),xycoords='axes fraction',fontsize=20)
ax1.annotate(r'$f_1(n_1,n_2)=0$',xy=(0.6,0.27),xycoords='axes fraction',fontsize=20)

ax1.set_xticklabels([])
ax1.set_yticklabels([])
plt.xlabel(r"$n_1$",fontsize=20)
plt.ylabel(r"$n_2$",fontsize=20)


#fig1, ax1 = plt.subplots()
#ax1.streamplot(X, Y, U, V,density=0.5)
#plt.xlim([0,2])
#plt.ylim([0,2])
#
#x=np.linspace(0,2,100)
#plt.plot(x,1-0.5*x,'k',linewidth=2)
#plt.plot(x,2-2*x,'k',linewidth=2)