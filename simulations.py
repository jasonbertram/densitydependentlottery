# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:01:54 2016

@author: jbertram
"""

from numpy import *
from pylab import *
import bisect

def expCisim(Ms,U,cs):
    scatter=[[0 for x in range(len(Ms))] for y in range(U)];
    for i in range(len(Ms)):
        for y in xrange(Ms[i]):
            cell=int(U*rand());
            scatter[cell][i]=scatter[cell][i]+1;
    
    winsnocomp=zeros(len(Ms)); wins1=zeros(len(Ms)); wins2=zeros(len(Ms));
    comp=zeros(U);
    for i in range(U):
        comp[i]=sum(scatter[i])
        if comp[i]>0:            
            lotterycmf=cumsum(array(scatter[i])*cs)
            victor=bisect.bisect(lotterycmf,rand()*lotterycmf[-1])
            
            if scatter[i][victor]==1 and comp[i]==1:
                winsnocomp[victor]=winsnocomp[victor]+1
            elif scatter[i][victor]==1:
                wins1[victor]=wins1[victor]+1
            else: 
                wins2[victor]=wins2[victor]+1
        
    return winsnocomp,wins1,wins2

def deltnu(Ms,U,cs,i):
    return Ms[i]*exp(-sum(Ms)/float(U))

def deltnr(Ms,U,cs,i):
    l=Ms[i]/float(U);notl=sum(Ms)/float(U)-l; L=float(l+notl);
    #return Ms[i]*exp(-l)*(1-exp(-notl))*cs[i]/(sum([Ms[j]*cs[j]/(U*(1-exp(-notl))) for j in range(len(Ms)) if (j!=i)])+cs[i])
    return Ms[i]*exp(-l)*(1-exp(-notl))*cs[i]/(sum([Ms[j]*cs[j]/(U*(1-(1+L)*exp(-L))) for j in range(len(Ms)) if (j!=i)])*(L-1+exp(-L))/notl+cs[i])
    
#def deltna(Ms,U,cs,i):
#    l=Ms[i]/float(U);
#    cbar=sum([Ms[j]*cs[j] for j in range(len(Ms))])/sum(Ms)
#    Da=cbar*(1-exp(-l))/(cbar*sum(Ms)/float(U)+cs[i]*l**2*exp(-l)/(1-(1+l)*exp(-l)))
#    return Ms[i]*Da*cs[i]/cbar

def deltna(Ms,U,cs,i):
    l=Ms[i]/float(U); L=sum(Ms)/float(U)
    lfac=(1-exp(-l))/(1-(1+l)*exp(-l))
    ljfac=(L*(1-exp(-L))/(1-(1+L)*exp(-L))-l*(1-exp(-l))/(1-(1+l)*exp(-l)))/(L-l)
    return Ms[i]*(1-exp(-l))*cs[i]/(lfac*l*cs[i]+ljfac*sum([Ms[j]*cs[j]/float(U) for j in range(len(Ms)) if (j!=i)]))
    
def deltnp(Ms,U,cs):
    return array([deltnu(Ms,U,cs,j)+deltnr(Ms,U,cs,j)+deltna(Ms,U,cs,j) for j in range(len(Ms))])

#Testing
#=================================================================
#Ms=[10000,100000,100000]; M=float(sum(Ms)); cs=array([1.,1.,5.]); G=len(Ms);
#print(expCisim(Ms,210000,cs))
#print(deltna(Ms,120000,cs,0))

#Simulation comparison figure
#================================================================
Ms=map(int,array([10000,90000]))
M=float(sum(Ms)); cs=array([1.5,1.]); G=len(Ms);
maxl=2; Us=array([int(Ms[0]/x) for x in linspace(0.001,maxl,1000)]);
exact=array([expCisim(Ms,U,cs) for U in Us]);

nocomp=array([deltnu(Ms,U,cs,0) for U in Us]);
comp1=array([deltnr(Ms,U,cs,0) for U in Us]);
comp2=array([deltna(Ms,U,cs,0) for U in Us]);

subplot(221)
plot(float(Ms[0])/Us,(exact[:,0,0]+exact[:,1,0]+exact[:,2,0])/Ms[0],'k.',markersize=1.,label="Simulation")
plot(float(Ms[0])/Us,(nocomp+comp1+comp2)/Ms[0],'k',label=r"$\Delta_+ n_1/m_1$")
#ylim([0,1000])
xlim([0,maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
gca().set_xticklabels(['0','','1','','2','','3'])
xlabel(r"$l_1$",fontsize=14)
ylabel(r"$\Delta_+ n_1/m_1$",fontsize=14)

gca().annotate(r'$(a)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(222)
plot(float(Ms[0])/Us,exact[:,0,0]/Ms[0],'k.',markersize=1, label="Simulation")
plot(float(Ms[0])/Us,exact[:,1,0]/Ms[0],'k.',markersize=1)
plot(float(Ms[0])/Us,exact[:,2,0]/Ms[0],'k.',markersize=1)

plot(float(Ms[0])/Us,nocomp/Ms[0],'k',label=r"$e^{-L}$")
plot(float(Ms[0])/Us,comp1/Ms[0],'r',label=r"$R_1 c_1/\overline{c}$")
plot(float(Ms[0])/Us,comp2/Ms[0],'b',label=r"$A_1 c_1/\overline{c}$")
#ylim([0,1000])
xlim([0,maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_1$",fontsize=14)
gca().set_xticklabels(['0','','1','','2','','3'])
gca().set_yticklabels([])
gca().annotate(r'$(b)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

nocomp=array([deltnu(Ms,U,cs,1) for U in Us]);
comp1=array([deltnr(Ms,U,cs,1) for U in Us]);
comp2=array([deltna(Ms,U,cs,1) for U in Us]);
    
subplot(223)
plot(float(Ms[1])/Us,(exact[:,0,1]+exact[:,1,1]+exact[:,2,1])/Ms[1],'k.',markersize=1,label="Simulation")
plot(float(Ms[1])/Us,(nocomp+comp1+comp2)/Ms[1],'k',label=r"$\Delta_+ n_2/m_2$")
#ylim([0,10000])
xlim([0,3*maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_2$",fontsize=14)
ylabel(r"$\Delta_+ n_2/m_2$",fontsize=14)
gca().annotate(r'$(c)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

subplot(224)
plot(float(Ms[1])/Us,exact[:,0,1]/Ms[1],'k.',markersize=1,label="Simulation")
plot(float(Ms[1])/Us,exact[:,1,1]/Ms[1],'k.',markersize=1)
plot(float(Ms[1])/Us,exact[:,2,1]/Ms[1],'k.',markersize=1)

plot(float(Ms[1])/Us,nocomp/Ms[1],'k',label=r"$e^{-L}$")
plot(float(Ms[1])/Us,comp1/Ms[1],'r',label=r"$R_2 c_2/\overline{c}$")
plot(float(Ms[1])/Us,comp2/Ms[1],'b',label=r"$A_2 c_2/\overline{c}$")
#ylim([0,10000])
xlim([0,3*maxl])
gca().xaxis.set_label_coords(0.5, -0.09)
xlabel(r"$l_2$",fontsize=14)
gca().set_yticklabels([])
gca().annotate(r'$(d)$',xy=(0.9,0.9),xycoords='axes fraction',fontsize=12)
legend(loc='upper center',prop={'size':10})

savefig('/home/jbertram/Dropbox/rcKderivation/simulationcomparison.pdf')


bs=array([0.6,1.05]);cs=array([20.,1.]);ds=array([0.5,0.5]);Ks=map(float,array([100000,100000]))

n0=array([20000,60000]); nt=[n0];
for i in range(10000):
    nt.append(nt[-1]+deltni(nt[-1],cs,Ks))

nt=array(nt)

print("first",bs[0]*cs[0]/ds[0])
print("second",bs[1]*cs[1]/ds[1])

plot(nt[:,0]+nt[:,1])
plot(nt[:,0],'k')
plot(nt[:,1],'k--')

#coexistence 
#=================================================================
K=100000; bs=array([0.4,.2]); cs=array([1.,1.]); ds=array([.2,.1]);
ns0=(1-ds)*array([1000,99000]);gens=100;

nhistapprox=[]; Us=[]; Rs=[]; As=[]; ns=ns0
for i in range(gens):
    Us.append(K-int(sum(ns)))
    #Rs.append([deltnr(bs*ns,Us[-1],cs,0),deltnr(bs*ns,Us[-1],cs,1)])
    #As.append([deltna(bs*ns,Us[-1],cs,0),deltna(bs*ns,Us[-1],cs,1)])
    delt=deltnp(bs*ns,Us[-1],cs)
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
#    temp=expCisim(map(int,bs*ns),Us[-1],cs)
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
