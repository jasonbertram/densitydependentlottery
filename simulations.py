# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:01:54 2016

@author: jbertram
"""

from numpy import *
from pylab import *
import bisect

def deltnplussim(m,c,U):
    scatter=[[0 for x in range(len(m))] for y in range(U)];
    for i in range(len(m)):
        for y in xrange(m[i]):
            cell=int(U*rand());
            scatter[cell][i]=scatter[cell][i]+1;
    
    winsnocomp=zeros(len(m)); wins1=zeros(len(m)); wins2=zeros(len(m));
    comp=zeros(U);
    for i in range(U):
        comp[i]=sum(scatter[i])
        if comp[i]>0:            
            lotterycmf=cumsum(array(scatter[i])*c)
            victor=bisect.bisect(lotterycmf,rand()*lotterycmf[-1])
            
            if scatter[i][victor]==1 and comp[i]==1:
                winsnocomp[victor]=winsnocomp[victor]+1
            elif scatter[i][victor]==1:
                wins1[victor]=wins1[victor]+1
            else: 
                wins2[victor]=wins2[victor]+1
        
    return winsnocomp,wins1,wins2

def R(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m); 
    return cbar*exp(-l)*(1-exp(-(L-l)))\
            /(c + (L-1+exp(-L))/(1-(1+L)*exp(-L))*(cbar*L-c*l)/(L-l))

def A(m,c,U):
    l=m/float(U)
    L=sum(l)
    cbar=sum(m*c)/sum(m)
    return (1-exp(-l))\
            /((1-exp(-l))/(1-(1+l)*exp(-l))*c*l\
            +(L*(1-exp(-L))/(1-(1+L)*exp(-L))-l*(1-exp(-l))/(1-(1+l)*exp(-l)))/(L-l)*(sum(c*l)-c*l))
    
def deltnplus(m,c,U):
    L=sum(m)/U
    cbar=sum(m*c)/sum(m)
    return m*(exp(-L)+(R(m,c,U)+A(m,c,U))*c/cbar)
    
def deltnplusclassic(m,c,U):
    return U*m*c/sum(m*c)
    

#Simulation comparison figure
#================================================================
m=array([10000,90000])
M=float(sum(m)) 
c=array([1.5,1.])
cbar=sum(m*c)/sum(m)
maxl=2; Us=array([int(m[0]/x) for x in linspace(0.001,maxl,100)]);
exact=array([deltnplussim(m,c,U) for U in Us]);
totalclassic=array([deltnplusclassic(m,c,U) for U in Us]);

#First genotype
nocomp=array([exp(-M/U) for U in Us]);
comp1=array([R(m,c,U)[0]*c[0]/cbar for U in Us]);
comp2=array([A(m,c,U)[0]*c[0]/cbar for U in Us]);

subplot(221)
plot(float(m[0])/Us,(exact[:,0,0]+exact[:,1,0]+exact[:,2,0])/m[0],'k.',markersize=1.,label="Simulation")
plot(float(m[0])/Us,(nocomp+comp1+comp2),'k',label=r"$\Delta_+ n_1/m_1$")
plot(float(m[0])/Us,totalclassic[:,0]/m[0],'k--',label=r"Classic lottery")
xlim([0,maxl])
ylim([0,1])
gca().xaxis.set_label_coords(0.5, -0.09)
gca().set_xticklabels(['0','','1','','2','','3'])
xlabel(r"$l_1$",fontsize=14)
ylabel(r"$\Delta_+ n_1/m_1$",fontsize=14)

gca().annotate(r'$(a)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(222)
plot(float(m[0])/Us,exact[:,0,0]/m[0],'k.',markersize=1, label="Simulation")
plot(float(m[0])/Us,exact[:,1,0]/m[0],'k.',markersize=1)
plot(float(m[0])/Us,exact[:,2,0]/m[0],'k.',markersize=1)

plot(float(m[0])/Us,nocomp,'k',label=r"$e^{-L}$")
plot(float(m[0])/Us,comp1,'r',label=r"$R_1 c_1/\overline{c}$")
plot(float(m[0])/Us,comp2,'b',label=r"$A_1 c_1/\overline{c}$")
#ylim([0,1000])
xlim([0,maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_1$",fontsize=14)
gca().set_xticklabels(['0','','1','','2','','3'])
gca().set_yticklabels([])
gca().annotate(r'$(b)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

#Second genotype
nocomp=array([exp(-M/U) for U in Us]);
comp1=array([R(m,c,U)[1]*c[1]/cbar for U in Us]);
comp2=array([A(m,c,U)[1]*c[1]/cbar for U in Us]);
    
subplot(223)
plot(float(m[1])/Us,(exact[:,0,1]+exact[:,1,1]+exact[:,2,1])/m[1],'k.',markersize=1,label="Simulation")
plot(float(m[1])/Us,(nocomp+comp1+comp2),'k',label=r"$\Delta_+ n_2/m_2$")
plot(float(m[1])/Us,totalclassic[:,1]/m[1],'k--',label=r"Classic lottery")
xlim([0,3*maxl])
ylim([0,1])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_2$",fontsize=14)
ylabel(r"$\Delta_+ n_2/m_2$",fontsize=14)
gca().annotate(r'$(c)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(224)
plot(float(m[1])/Us,exact[:,0,1]/m[1],'k.',markersize=1,label="Simulation")
plot(float(m[1])/Us,exact[:,1,1]/m[1],'k.',markersize=1)
plot(float(m[1])/Us,exact[:,2,1]/m[1],'k.',markersize=1)

plot(float(m[1])/Us,nocomp,'k',label=r"$e^{-L}$")
plot(float(m[1])/Us,comp1,'r',label=r"$R_2 c_2/\overline{c}$")
plot(float(m[1])/Us,comp2,'b',label=r"$A_2 c_2/\overline{c}$")
#ylim([0,10000])
xlim([0,3*maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_2$",fontsize=14)
gca().set_yticklabels([])
gca().annotate(r'$(d)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

savefig('/home/jbertram/repos/densitydependentlottery/simulationcomparison.pdf')

#Seasonal population fluctuations
#=================================================================

def b(t):
    return bparam[0]+bparam[1]*(1+sin(2*pi*t/g))
    
def d(t):
    return dparam[0]+dparam[1]*(1+sin(2*pi*t/g))

n=array([5000,5000])
c=array([1.,1.])
T=100000
g=2.     #generations per seasonal cycle
bparam=array([[0.1,0.1],[0.1,0.1]])
dparam=array(array([[0.2,0.2],[0.,0.]]))
    
nhist=[]
totaltime=20
for t in range(totaltime):
    print T-sum(n)
    n=n+deltnplus(b(t)*n,c,T-sum(n))
    nhist.append(list(n))

#coexistence 
#=================================================================
K=100000; bs=array([0.4,.2]); c=array([1.,1.]); ds=array([.2,.1]);
ns0=(1-ds)*array([1000,99000]);gens=100;

nhistapprox=[]; Us=[]; Rs=[]; As=[]; ns=ns0
for i in range(gens):
    Us.append(K-int(sum(ns)))
    #Rs.append([deltnr(bs*ns,Us[-1],c,0),deltnr(bs*ns,Us[-1],c,1)])
    #As.append([deltna(bs*ns,Us[-1],c,0),deltna(bs*ns,Us[-1],c,1)])
    delt=deltnplus(bs*ns,Us[-1],c)
    ns=ns+delt-ns*ds;
    nhistapprox.append(ns)

nhistapprox=array(nhistapprox)
#Rs=array(Rs)
#As=array(As)

figure()
plot(nhistapprox[:,0])
plot(nhistapprox[:,1])

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
#nhist=array(nhist)
#Rsexact=array(Rsexact)
#Asexact=array(Asexact)
#
#plot(nhist[:,0])
#plot(nhist[:,1])

#figure()
#plot(Rs[:,0],'k')
#plot(Rs[:,1],'b')
#plot(As[:,0],'k--')
#plot(As[:,1],'b--')
#
#plot(Rsexact[:,0],'k.')
#plot(Rsexact[:,1],'b.')
#plot(Asexact[:,0],'kx')
#plot(Asexact[:,1],'bx')



figure(figsize=[4,3])

fill_between(range(gens),0,nhistapprox[:,0]/(nhistapprox[:,0]+nhistapprox[:,1]))
ylim([0,1])

#gca().set_xticklabels([]);
#gca().set_yticklabels([]);
gca().xaxis.set_label_coords(0.5, -0.075)
xlabel('Generations')
gca().yaxis.set_label_coords(-0.1, 0.5)
ylabel('Frequency')

savefig('/home/jbertram/Dropbox/rcKderivation/coex.pdf')
