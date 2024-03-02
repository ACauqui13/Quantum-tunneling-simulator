# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 14:58:00 2023
assymetric biased double barrier
@author: Alvaro Cauqui Diaz

"""
import numpy as np
import math
import matplotlib.pyplot as plt
i=1j
me=9.109e-31
m=0.067*me
def V3(x_):
    x= x_ * 1e-9
    if x<1e-9:
        return 0
    elif 1e-9<=x<5e-9:
        m1 = -5e7
        b1=0.55
        return m1*x +b1
    elif 5e-9 <= x < 9e-9:
        m2 = -5e7
        b2=0.05
        return m2*x + b2
    elif 9e-9<=x<12e-9:
        m3=-5e7
        b3=1.15
        return m3*x +b3
    else:
        return -0.63
    
x0=[x for x in np.arange(0,13,0.001)]
Ve=[V3(f) for f in x0]
plt.plot(x0,Ve)
plt.hlines(y=0.4, xmin=-1, xmax=13, color='black', linestyles='--')
plt.hlines(y=0.6, xmin=-1, xmax=13, color='black', linestyles='--')
plt.xlabel('x(nm)')
plt.ylabel('V(eV)')

#%%
def evtoj(ev):
    j=ev*1.602e-19
    
    return j
def nmtom(nm):
    _=nm*10**-9
    return _
def K(x2,x4,k1,k3,k5,w2,w3,w4):
    
    phi1=k3*w3
    phi2=np.arctan(x2/k1)
    phi3=np.arctan(x2/k3)
    phi4=np.arctan(x4/k3)
    phi5=np.arctan(x4/k5)
    return np.exp((x2*w2)+(x4*w4))*(np.exp(i*(-phi1+phi2+phi3+phi4+phi5))-np.exp(i*(phi1+phi2-phi3-phi4+phi5)))\
        +np.exp((x2*w2)-(x4*w4))*(-np.exp(i*(-phi1+phi2+phi3-phi4-phi5))+np.exp(i*(phi1+phi2-phi3+phi4-phi5))) \
            +np.exp((-x2*w2)+(x4*w4))*(-np.exp(i*(-phi1-phi2-phi3+phi4+phi5))+np.exp(i*(phi1-phi2+phi3-phi4+phi5))) \
                +np.exp((-x2*w2)+(-x4*w4))*(np.exp(i*(-phi1-phi2-phi3-phi4-phi5))-np.exp(i*(phi1-phi2+phi3+phi4-phi5)))
def k(Ei,Vx,m=m):
    hbar=1.054571817e-34
    #Vx=V(x)*1.6e-19
    return np.sqrt((2*m*(Ei-Vx)))/hbar
def x(Ei,Vx,m=m):
    hbar=1.054571817e-34
    #Vx=V(x)*1.6e-19
    return np.sqrt((2*m*(Vx-Ei)))/hbar
def T(x2,x4,k1,k3,k5,w2,w3,w4):
    K_=K(x2,x4,k1,k3,k5,w2,w3,w4)
    Ksq=K_.conjugate()*K_
    return (((2**8)*k1*(x2**2)*(k3**2)*(x4**2)*k5)) \
        /((Ksq.real)*(k1**2+x2**2)*(x2**2+k3**2)*(k3**2+x4**2)*(x4**2+k5**2))

def Tc(E):
    Ej=evtoj(E)
    k1=k(Ej,0,me)
    k3=k(Ej,(-0.3)*1.6e-19,m)
    k5=k(Ej,((-0.67)*1.6e-19),m)
    x4=x(Ej,(0.6)*1.6e-19,m)
    x2=x(Ej,(0.4)*1.6e-19,m)
    w2=4e-9
    w3=4e-9
    w4=3e-9
    return T(x2,x4,k1,k3,k5,w2,w3,w4)

Ee=[i for i in np.arange(0.001,0.399,0.001)]

Tg=[Tc(g) for g in Ee]
plt.figure()
plt.plot(Ee,Tg,color='black')
plt.xlabel('E(eV)')
plt.ylabel('T')
    
    