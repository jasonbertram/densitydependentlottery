# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 12:26:40 2017

@author: jbertram
"""

#(n1,n2) vector field
#=================================================================
T=100000
x,y=np.linspace(1e-5,1.,401),np.linspace(1e-5,1.,401)
xgrid, ygrid=np.meshgrid(x,y)
X,Y=np.zeros([len(x),len(y)]),np.zeros([len(x),len(y)])
b=np.array([5.,5.])
c=np.array([5.,1.])
d=np.array([1.1,1.1])

for i in xrange(len(x)):
    for j in xrange(len(y)):
        n=T*np.array([x[i],y[j]])
        U=T-sum(n)
        X[j,i],Y[j,i]=((1+deltnplus(b*n*U/float(T),c,U)/n)/d-1)*n


fig1, (ax2, ax1) = plt.subplots(ncols=2,figsize=[8,4])

f1null=np.array([[x[[_ for _ in range(len(x)) if X[i,_]>0][-1]],y[i]] for i in range(1,len(y)) if max(X[i])>0])
f2null=np.array([[x[[_ for _ in range(len(x)) if Y[i,_]>0][-1]],y[i]] for i in range(1,len(y)) if max(Y[i])>0])
ax1.plot(f1null[:,0],f1null[:,1],'k--',linewidth=2.)
ax1.plot(f2null[:,0],f2null[:,1],'k--',linewidth=2.)
ax1.plot([f2null[-1,0],f1null[0,0]],[f2null[-1,1],f1null[0,1]],'k',linewidth=2,alpha=0.8)
ax1.set_xlim([0,f1null[0,0]])
ax1.set_ylim([0,f1null[-1,1]])

ax1.annotate(r'$(b)$',xy=(0.9,0.93),xycoords='axes fraction',fontsize=16)

seedpoints=[[x[_]/2,y[200-_]/2] for _ in range(0,201,10)]
seedpoints=np.concatenate([[[_,0.2-_] for _ in np.arange(0.02,0.2,0.02)],\
            [[float(_)/10,1.-_*0.1] for _ in range(0,11) if _ not in []]])

ax1.streamplot(x, y, X, Y,start_points=seedpoints,linewidth=1.,density=100)

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_xlabel(r"$n_1$",fontsize=20)
ax1.set_ylabel(r"$n_2$",fontsize=20)

x,y=np.linspace(0,1.5,1000),np.linspace(0,1.5,1000)
X, Y = np.meshgrid(x,y)
alphas=np.array([[1.,.9],[1.2,1.]])
U = 1*(1-alphas[0,0]*X-alphas[0,1]*Y)*X
V = 1*(1-alphas[1,1]*Y-alphas[1,0]*X)*Y

z1=np.linspace(0,1,10)
ax2.plot(z1,(1-alphas[0,0]*z1)/alphas[0,1],'k--',linewidth=2,label="Nullclines")
ax2.plot(z1,(1-alphas[1,0]*z1)/alphas[1,1],'k--',linewidth=2)
ax2.plot([0,1/alphas[0,0]],[1/alphas[1,1],0],'k',linewidth=2,alpha=0.8,label=r"$N=K_{11}=K_{22}$")
ax2.set_xlim([0,1/alphas[0,0]])
ax2.set_ylim([0,1/alphas[0,1]])

ax2.legend(loc='upper center',prop={'size':11})

seedpoints=np.concatenate([[[_,0.2-_] for _ in np.arange(0.02,0.2,0.02)],\
            [[float(_)/10,1.-_*0.1] for _ in range(0,11) if _ not in []]])
ax2.streamplot(x, y, U, V,density=100,start_points=seedpoints,linewidth=1.)

#ax2.annotate(r'$\frac{dn_2}{dt}=0$',xy=(0.15,0.81),xycoords='axes fraction',fontsize=16)
#ax2.annotate(r'$\frac{dn_1}{dt}=0$',xy=(0.62,0.28),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$K_{21}$',xy=(0.82,-0.07),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$K_{11}$',xy=(0.95,-0.07),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$K_{22}$',xy=(-0.13,0.89),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$K_{12}$',xy=(-0.13,.97),xycoords='axes fraction',fontsize=16)
ax2.annotate(r'$(a)$',xy=(0.9,0.93),xycoords='axes fraction',fontsize=16)

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xlabel(r"$n_1$",fontsize=20)
ax2.set_ylabel(r"$n_2$",fontsize=20)

plt.tight_layout()

plt.savefig('/home/jbertram/repos/densitydependentlottery/LVvslottery.pdf',bbox="tight")