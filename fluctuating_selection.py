# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:01:54 2016

@author: jbertram
"""


#Seasonal population fluctuations
#=================================================================

def b(t):
    #return bparam[0]+bparam[1]*(1+sin(2*pi*t/g))
    if t%g<g/2:
        return bparam[0]
    else:
        return bparam[1]
    
def d(t):
    #return dparam[0]+dparam[1]*(1+cos(2*pi*t/g))
    if t%g<g/2:
        return dparam[1]
    else:
        return dparam[0]


c=np.array([1.,1.])
T=100000
g=40.     #iterations per seasonal cycle
bparam=np.array([[0.0,0.0],[2.,.5]])
dparam=np.array([[0.0,0.0],[0.2,0.1]])
totaltime=100

n=np.array([60000.,25000.])
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n = n + deltnplus((b(t)*n*U/T).astype(np.int),c,U.astype(np.int))-d(t)*n
    nhist.append(list(n))
    
fig=plt.subplot(331)
for t in range(0,100,40):
    fig.plot(range(t,t+20),[dparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
    fig.plot(range(t,t+20),[dparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
    
fig.fill([20,40,40,20],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.fill([60,40+40,40+40,60],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.plot(range(20,40),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
fig.plot(range(20,40),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
fig.plot(range(60,80),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2,label=r"$b$ spec.")
fig.plot(range(60,80),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2,label=r"$d$ spec.")

fig.set_xticklabels([])
fig.set_xlim([0,100])
fig.set_ylim([0,0.7])
plt.ylabel(r"$b(t), d(t)$",fontsize=14)
fig.annotate(r'$(a)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)
fig.annotate(r'$d$',xy=(0.09,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.29,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.49,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.69,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.89,0.4),xycoords='axes fraction',fontsize=16)

fig=plt.subplot(334)
plt.stackplot(range(totaltime+1),np.array(nhist)[:,1],np.array(nhist)[:,0],colors=['b','k'])
fig.set_xticklabels([])
fig.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
plt.ylabel(r"Abundance$/T$",fontsize=14)
fig.annotate(r'$(b)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)

fig=plt.subplot(337)
plt.stackplot(range(totaltime+1),np.array(nhist)[:,1]/np.sum(nhist,1),np.array(nhist)[:,0]/np.sum(nhist,1),colors=['b','k'])
plt.ylim([0,1])
plt.xlabel(r"$t$",fontsize=14)
plt.ylabel("Frequency",fontsize=14)
fig.annotate(r'$(c)$',xy=(0.01,0.9),xycoords='axes fraction',color='white',fontsize=12)

#Classic lottery same
n=np.array([60000.,25000.])
nhistclassic=[n];
for t in range(totaltime):
    U=T-sum(n)
    n = n + deltnplusclassicboundedm((b(t)*n*U/T).astype(np.int),c,U.astype(np.int))-d(t)*n
    nhistclassic.append(list(n))

fig=plt.subplot(332)
for t in range(0,100,40):
    fig.plot(range(t,t+20),[dparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
    fig.plot(range(t,t+20),[dparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
    
fig.fill([20,40,40,20],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.fill([60,40+40,40+40,60],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.plot(range(20,40),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
fig.plot(range(20,40),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
fig.plot(range(60,80),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2,label=r"$b$ spec.")
fig.plot(range(60,80),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2,label=r"$d$ spec.")

fig.set_xticklabels([])
fig.set_xlim([0,100])
fig.set_ylim([0,0.7])
fig.annotate(r'$(d)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)
fig.annotate(r'$d$',xy=(0.09,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.29,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.49,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.69,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.89,0.4),xycoords='axes fraction',fontsize=16)

fig=plt.subplot(335)
plt.stackplot(range(totaltime+1),np.array(nhistclassic)[:,1],np.array(nhistclassic)[:,0],colors=['b','k'])
fig.set_xticklabels([])
fig.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
fig.annotate(r'$(e)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)

fig=plt.subplot(338)
plt.stackplot(range(totaltime+1),np.array(nhistclassic)[:,1]/np.sum(nhistclassic,1),np.array(nhistclassic)[:,0]/np.sum(nhistclassic,1),colors=['b','k'])
plt.ylim([0,1])
plt.xlabel(r"$t$",fontsize=14)
fig.annotate(r'$(f)$',xy=(0.01,0.9),xycoords='axes fraction',color='white',fontsize=12)

#Classic lottery coex
c=np.array([1.,1.])
T=100000
g=40.     #iterations per seasonal cycle
bparam=np.array([[0.0,0.0],[.5,.0421]])
dparam=np.array([[0.0,0.0],[0.2,0.1]])
totaltime=100

n=np.array([60000.,25000.])
nhistclassic=[n];
for t in range(totaltime):
    U=T-sum(n)
    n = n + deltnplusclassicboundedm((b(t)*n*U/T).astype(np.int),c,U.astype(np.int))-d(t)*n
    nhistclassic.append(list(n))
    
fig=plt.subplot(333)
for t in range(0,100,40):
    fig.plot(range(t,t+20),[dparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
    fig.plot(range(t,t+20),[dparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
    
fig.fill([20,40,40,20],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.fill([60,40+40,40+40,60],[0,0,1,1],'g',alpha=0.3,edgecolor='none')
fig.plot(range(20,40),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2)
fig.plot(range(20,40),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2)
fig.plot(range(60,80),[bparam[1][0] for _ in range(t,t+20)],'k.',lw=2,label=r"$b$ specialist")
fig.plot(range(60,80),[bparam[1][1] for _ in range(t,t+20)],'b.',lw=2,label=r"$d$ specialist")

fig.set_xticklabels([])
fig.set_yticklabels([])
fig.set_xlim([0,100])
fig.set_ylim([0,0.7])
fig.annotate(r'$(g)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)
fig.annotate(r'$d$',xy=(0.09,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.29,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.49,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.69,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.89,0.4),xycoords='axes fraction',fontsize=16)
plt.legend(loc='upper right',prop={'size':10})

fig=plt.subplot(336)
plt.stackplot(range(totaltime+1),np.array(nhistclassic)[:,1],np.array(nhistclassic)[:,0],colors=['b','k'])
fig.set_xticklabels([])
fig.set_yticklabels([])
fig.annotate(r'$(h)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)

fig=plt.subplot(339)
plt.stackplot(range(totaltime+1),np.array(nhistclassic)[:,1]/np.sum(nhistclassic,1),np.array(nhistclassic)[:,0]/np.sum(nhistclassic,1),colors=['b','k'])
plt.ylim([0,1])
plt.xlabel(r"$t$",fontsize=14)
fig.set_yticklabels([])
fig.annotate(r'$(i)$',xy=(0.01,0.9),xycoords='axes fraction',color='white',fontsize=12)

plt.savefig('/home/jbertram/repos/densitydependentlottery/fluctuatingselection.pdf')

    
#n=np.array([6000.,2500.])
#nsimhist=[n]
#for t in range(totaltime):
#    U=T-sum(n)
#    n = n + sum(deltnplussim((b(t)*n*U/T).astype(np.int),c,U.astype(np.int)))-d(t)*n
#    nsimhist.append(list(n))

#coexistence 
##=================================================================
#K=100000; bs=np.array([0.4,.2]); c=np.array([1.,1.]); ds=np.array([.2,.1]);
#ns0=(1-ds)*np.array([1000,99000]);gens=100;
#
#nhistapprox=[]; Us=[]; Rs=[]; As=[]; ns=ns0
#for i in range(gens):
#    Us.append(K-int(sum(ns)))
#    #Rs.append([deltnr(bs*ns,Us[-1],c,0),deltnr(bs*ns,Us[-1],c,1)])
#    #As.append([deltna(bs*ns,Us[-1],c,0),deltna(bs*ns,Us[-1],c,1)])
#    delt=deltnplus(bs*ns,Us[-1],c)
#    ns=ns+delt-ns*ds;
#    nhistapprox.append(ns)
#
#nhistapprox=np.array(nhistapprox)
##Rs=np.array(Rs)
##As=np.array(As)
#
#plt.figure()
#plt.plot(nhistapprox[:,0])
#plt.plot(nhistapprox[:,1])

#nhist=[]; Us=[]; ns=ns0; Rsexact=[];Asexact=[]
#for i in range(gens):
#    Us.append(K-int(sum(ns)))    
#    temp=deltnsim(map(int,bs*ns),Us[-1],c)
#    Rsexact.append(temp[1])
#    Asexact.append(temp[2])
#    delt=sum(temp,0)
#    ns=ns+delt-ns*ds;
#    nhist.append(ns)
#
#nhist=np.array(nhist)
#Rsexact=np.array(Rsexact)
#Asexact=np.array(Asexact)
#
#plt.plot(nhist[:,0])
#plt.plot(nhist[:,1])

#plt.figure()
#plt.plot(Rs[:,0],'k')
#plt.plot(Rs[:,1],'b')
#plt.plot(As[:,0],'k--')
#plt.plot(As[:,1],'b--')
#
#plt.plot(Rsexact[:,0],'k.')
#plt.plot(Rsexact[:,1],'b.')
#plt.plot(Asexact[:,0],'kx')
#plt.plot(Asexact[:,1],'bx')



#plt.figure(figsize=[4,3])
#
#fill_between(range(gens),0,nhistapprox[:,0]/(nhistapprox[:,0]+nhistapprox[:,1]))
#ylim([0,1])
#
##gca().set_xticklabels([]);
##gca().set_yticklabels([]);
#gca().xaxis.set_label_coords(0.5, -0.075)
#xlabel('Generations')
#gca().yaxis.set_label_coords(-0.1, 0.5)
#ylabel('Frequency')
#
#savefig('/home/jbertram/Dropbox/rcKderivation/coex.pdf')
