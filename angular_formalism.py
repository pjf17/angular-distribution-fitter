from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from sympy.physics.quantum.cg import CG
from sympy import S
import sympy as sp
import math
import numpy as np
from scipy.special import legendre
from scipy.special import lpmv

#================================================================================
#angular momementum coupling
#================================================================================
def kappa(k, l, lp):
    #modifier coefficients for linear polarization distribution
    if (k > l + lp or k == 0):
        return 0
    prefactor = -1*np.sqrt(np.math.factorial(k-2)/np.math.factorial(k+2))
    cgnum = CG(l,1,lp,1,k,2)
    cgden = CG(l,1,lp,-1,k,0)
    return prefactor*float(cgnum.doit())/float(cgden.doit())

def F_k(k, l, lp, If, Ii):
    #ordinary F-coefficient
    prefactor = (-1)**(If+Ii+1) * (2*k+1)**(1/2) * (2*l+1)**(1/2) * (2*lp+1)**(1/2) * (2*Ii+1)**(1/2)
    return prefactor * wigner_3j(l,lp,k,1,-1,0) * wigner_6j(l,lp,k,Ii,Ii,If)

def A_k(k, delta, l, lp, If, Ii):
    #A coeffients for the standard angular distribution formula
    return 1.0/(1.0 + delta**2) * (F_k(k,l,l,If,Ii) + 2*delta*F_k(k,l,lp,If,Ii) +  delta**2*F_k(k,lp,lp,If,Ii))

def Ap_k(k, delta, l, lp, If, Ii):
    #A' coeffients for the linear polarization distribution formula
    return 1.0/(1.0 + delta**2) * (-1.0*kappa(k,l,l)*F_k(k,l,l,If,Ii) + 2*delta*kappa(k,l,lp)*F_k(k,l,lp,If,Ii) +  delta**2*kappa(k,lp,lp)*F_k(k,lp,lp,If,Ii))

#================================================================================
#alignment
#================================================================================
def u_k(k, l, Ii, If):
    #coefficient for propagating alignment across multiple transitions !(w/o mixing)!
    prefactor = (-1)**(Ii+If+l+k) * ((2*Ii+1)*(2*If+1))**(1/2)
    return prefactor * wigner_6j(Ii,Ii,k,If,If,l)

def U_k(k, delta, l, lp, Ii, If):
    #coefficient for propagating alignment across multiple transitions !(w/ mixing)!
    return (u_k(k,l,Ii,If) + delta**2*u_k(k,lp,Ii,If))/(1+delta**2)

def substates(I):
    # creates a list of all substates for a given angular momentum
    Ms = []
    mm = -I
    for i in range(int(2*I)+1):
        Ms.append(mm)
        mm += 1
    
    return Ms

def calcPm(I,sigma,algntype):
    # calculates the substate population distribution under differen alignment type assumptions
    Pm = {}
    norm = 0
    subs = substates(I)
    width = sigma*I
    algnmodel = None
    if algntype: #prolate
        algnmodel = lambda x: sp.exp(-(I - abs(x))**2/(2*width**2))
    else: #oblate
        algnmodel = lambda x: sp.exp(-x**2/(2*width**2))

    if (sigma > 0):
        for m in subs:
            g = algnmodel(m)
            norm += g
            Pm.update({m: g})

        for m in subs:
            Pm[m] = Pm[m]/norm
    
    else: 
        #width is zero, so fully aligned, four (effectively 3) possible cases 
        #for prolate/oblate alignments and integer/half integer spins
        
        if algntype: #prolate is same for both half integer and integer
            for m in subs:
                if (m == abs(I)):
                    Pm.update({m: 0.5})
                else:
                    Pm.update({m: 0.0})
        
        else: #oblate
            if ((I - np.floor(I)) == 0): #integer spin
                for m in subs:
                    if (m != 0):
                        Pm.update({m: 0.0})
                    else:
                        Pm.update({m: 1.0})

            else: #half integer spin
                for m in subs:
                    if ((abs(m) != S(1)/2 )):
                        Pm.update({m: 0.0})
                    else:
                        Pm.update({m: 0.5})
    
    return Pm  

def B_k(k, Ii, sigma, algntype, Pm=None):
    #the coefficient for alignment 
    prefactor = (2*k+1)**(1/2) * (2*Ii + 1)**(1/2)
    sum = 0
    Ms = substates(Ii)

    if (Pm == None):
        Pm = calcPm(Ii,sigma,algntype)

    for m in Ms:
        sum += (-1)**(Ii + m)*wigner_3j(Ii,Ii,k,-m,m,0)*Pm[m]
    
    return prefactor*sum