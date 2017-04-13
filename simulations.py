# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:01:54 2016

@author: jbertram
"""

import numpy as np
import matplotlib.pyplot as plt
import bisect

def deltnplussim(m,c,U):
    scatter=[[0 for x in range(len(m))] for y in range(U)];
    for i in range(len(m)):
        for y in xrange(m[i]):
            cell=int(U*np.random.rand());
            scatter[cell][i]=scatter[cell][i]+1;
    
    winsnocomp=np.zeros(len(m)); wins1=np.zeros(len(m)); wins2=np.zeros(len(m));
    comp=np.zeros(U);
    for i in range(U):
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
        
    return np.array([winsnocomp,wins1,wins2])

def R(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m); 
    return cbar*np.exp(-l)*(1-np.exp(-(L-l)))\
            /(c + (L-1+np.exp(-L))/(1-(1+L)*np.exp(-L))*(cbar*L-c*l)/(L-l))

def A(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m)
    return (1-np.exp(-l))\
            /((1-np.exp(-l))/(1-(1+l)*np.exp(-l))*c*l\
            +(L*(1-np.exp(-L))/(1-(1+L)*np.exp(-L))-l*(1-np.exp(-l))/(1-(1+l)*np.exp(-l)))/(L-l)*(sum(c*l)-c*l))
    
def deltnplus(m,c,U):
    if sum(m)>0 and U>0:
        L=sum(m)/float(U)
        cbar=sum(m*c)/sum(m)
        return m*(np.exp(-L)+(R(m,c,U)+A(m,c,U))*c/cbar)
    else:
        return np.zeros(len(m))
    
def deltnplusclassic(m,c,U):
    if sum(m)>0:
        return U*m*c/sum(m*c)
    else:
        return np.zeros(len(m))

#Simulation comparison figure
#================================================================
m=np.array([10000,90000])
M=float(sum(m)) 
c=np.array([1.5,1.])
cbar=sum(m*c)/sum(m)
maxl=2; Us=np.array([int(m[0]/x) for x in linspace(0.001,maxl,100)]);
exact=np.array([deltnplussim(m,c,U) for U in Us]);
totalclassic=np.array([deltnplusclassic(m,c,U) for U in Us]);

#First genotype
nocomp=np.array([np.exp(-M/U) for U in Us]);
comp1=np.array([R(m,c,U)[0]*c[0]/cbar for U in Us]);
comp2=np.array([A(m,c,U)[0]*c[0]/cbar for U in Us]);

subplot(221)
plt.plot(float(m[0])/Us,(exact[:,0,0]+exact[:,1,0]+exact[:,2,0])/m[0],'k.',markersize=1.,label="Simulation")
plt.plot(float(m[0])/Us,(nocomp+comp1+comp2),'k',label=r"$\Delta_+ n_1/m_1$")
plt.plot(float(m[0])/Us,totalclassic[:,0]/m[0],'k--',label=r"Classic lottery")
xlim([0,maxl])
ylim([0,1])
gca().xaxis.set_label_coords(0.5, -0.09)
gca().set_xticklabels(['0','','1','','2','','3'])
xlabel(r"$l_1$",fontsize=14)
ylabel(r"$\Delta_+ n_1/m_1$",fontsize=14)

gca().annotate(r'$(a)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(222)
plt.plot(float(m[0])/Us,exact[:,0,0]/m[0],'k.',markersize=1, label="Simulation")
plt.plot(float(m[0])/Us,exact[:,1,0]/m[0],'k.',markersize=1)
plt.plot(float(m[0])/Us,exact[:,2,0]/m[0],'k.',markersize=1)

plt.plot(float(m[0])/Us,nocomp,'k',label=r"$e^{-L}$")
plt.plot(float(m[0])/Us,comp1,'r',label=r"$R_1 c_1/\overline{c}$")
plt.plot(float(m[0])/Us,comp2,'b',label=r"$A_1 c_1/\overline{c}$")
#ylim([0,1000])
xlim([0,maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_1$",fontsize=14)
gca().set_xticklabels(['0','','1','','2','','3'])
gca().set_yticklabels([])
gca().annotate(r'$(b)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

#Second genotype
nocomp=np.array([np.exp(-M/U) for U in Us]);
comp1=np.array([R(m,c,U)[1]*c[1]/cbar for U in Us]);
comp2=np.array([A(m,c,U)[1]*c[1]/cbar for U in Us]);
    
subplot(223)
plt.plot(float(m[1])/Us,(exact[:,0,1]+exact[:,1,1]+exact[:,2,1])/m[1],'k.',markersize=1,label="Simulation")
plt.plot(float(m[1])/Us,(nocomp+comp1+comp2),'k',label=r"$\Delta_+ n_2/m_2$")
plt.plot(float(m[1])/Us,totalclassic[:,1]/m[1],'k--',label=r"Classic lottery")
xlim([0,3*maxl])
ylim([0,1])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_2$",fontsize=14)
ylabel(r"$\Delta_+ n_2/m_2$",fontsize=14)
gca().annotate(r'$(c)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(224)
plt.plot(float(m[1])/Us,exact[:,0,1]/m[1],'k.',markersize=1,label="Simulation")
plt.plot(float(m[1])/Us,exact[:,1,1]/m[1],'k.',markersize=1)
plt.plot(float(m[1])/Us,exact[:,2,1]/m[1],'k.',markersize=1)

plt.plot(float(m[1])/Us,nocomp,'k',label=r"$e^{-L}$")
plt.plot(float(m[1])/Us,comp1,'r',label=r"$R_2 c_2/\overline{c}$")
plt.plot(float(m[1])/Us,comp2,'b',label=r"$A_2 c_2/\overline{c}$")
#plt.ylim([0,10000])
plt.xlim([0,3*maxl])
plt.gca().xaxis.set_label_coords(0.5, -0.09)
plt.xlabel(r"$l_2$",fontsize=14)
plt.gca().set_yticklabels([])
plt.gca().annotate(r'$(d)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
plt.legend(loc='upper center',prop={'size':10})

plt.savefig('/home/jbertram/repos/densitydependentlottery/simulationcomparison.pdf')

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
bparam=np.array([[0.0,0.0],[.5,.2]])
dparam=np.array([[0.0,0.0],[0.2,0.1]])
totaltime=100

n=np.array([60000.,25000.])
nhist=[n];
for t in range(totaltime):
    U=T-sum(n)
    n = n + deltnplus((b(t)*n*U/T).astype(np.int),c,U.astype(np.int))-d(t)*n
    nhist.append(list(n))
    
#n=np.array([60000.,25000.])
#nhistclassic=[n];
#for t in range(totaltime):
#    U=T-sum(n)
#    n = n + deltnplusclassic((b(t)*n*U/T).astype(np.int),c,U.astype(np.int))-d(t)*n
#    nhistclassic.append(list(n))
#
#plt.stackplot(range(totaltime+1),np.array(nhistclassic)[:,1],np.array(nhistclassic)[:,0],colors=['b','k'])
#plt.plot(np.array(nhistclassic)[:,0]/np.sum(nhistclassic,1),'k',lw=2)
    
    
#n=np.array([6000.,2500.])
#nsimhist=[n]
#for t in range(totaltime):
#    U=T-sum(n)
#    n = n + sum(deltnplussim((b(t)*n*U/T).astype(np.int),c,U.astype(np.int)))-d(t)*n
#    nsimhist.append(list(n))

fig=plt.subplot(311)
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
fig.set_xlim([0,100])
fig.set_ylim([0,0.7])
plt.ylabel(r"$b(t), d(t)$",fontsize=14)
fig.annotate(r'$(a)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)
fig.annotate(r'$d$',xy=(0.09,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.29,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.49,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$b$',xy=(0.69,0.4),xycoords='axes fraction',fontsize=16)
fig.annotate(r'$d$',xy=(0.89,0.4),xycoords='axes fraction',fontsize=16)
plt.legend(loc='upper center',prop={'size':10})

fig=plt.subplot(312)
plt.stackplot(range(totaltime+1),np.array(nhist)[:,1],np.array(nhist)[:,0],colors=['b','k'])
fig.set_xticklabels([])
fig.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
plt.ylabel(r"Abundance$/T$",fontsize=14)
fig.annotate(r'$(b)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)

fig=plt.subplot(313)
plt.stackplot(range(totaltime+1),np.array(nhist)[:,1]/np.sum(nhist,1),np.array(nhist)[:,0]/np.sum(nhist,1),colors=['b','k'])
plt.ylim([0,1])
plt.xlabel(r"$t$",fontsize=14)
plt.ylabel("Frequency",fontsize=14)
fig.annotate(r'$(c)$',xy=(0.01,0.9),xycoords='axes fraction',fontsize=12)

plt.savefig('/home/jbertram/repos/densitydependentlottery/fluctuatingselection.pdf')


#coexistence 
#=================================================================
K=100000; bs=np.array([0.4,.2]); c=np.array([1.,1.]); ds=np.array([.2,.1]);
ns0=(1-ds)*np.array([1000,99000]);gens=100;

nhistapprox=[]; Us=[]; Rs=[]; As=[]; ns=ns0
for i in range(gens):
    Us.append(K-int(sum(ns)))
    #Rs.append([deltnr(bs*ns,Us[-1],c,0),deltnr(bs*ns,Us[-1],c,1)])
    #As.append([deltna(bs*ns,Us[-1],c,0),deltna(bs*ns,Us[-1],c,1)])
    delt=deltnplus(bs*ns,Us[-1],c)
    ns=ns+delt-ns*ds;
    nhistapprox.append(ns)

nhistapprox=np.array(nhistapprox)
#Rs=np.array(Rs)
#As=np.array(As)

plt.figure()
plt.plot(nhistapprox[:,0])
plt.plot(nhistapprox[:,1])

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



plt.figure(figsize=[4,3])

fill_between(range(gens),0,nhistapprox[:,0]/(nhistapprox[:,0]+nhistapprox[:,1]))
ylim([0,1])

#gca().set_xticklabels([]);
#gca().set_yticklabels([]);
gca().xaxis.set_label_coords(0.5, -0.075)
xlabel('Generations')
gca().yaxis.set_label_coords(-0.1, 0.5)
ylabel('Frequency')

savefig('/home/jbertram/Dropbox/rcKderivation/coex.pdf')
