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
                if (abs(m) == I):
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
                    if ((abs(m) != S(1/2) )):
                        Pm.update({m: 0.0})
                    else:
                        Pm.update({m: 0.5})
    
    return Pm  

def B_k(k, **kwargs):
    #the coefficient for alignment
    #option 1: quantization axis is established by the beam direction 
    #PARAMETERS: either "sigma" and "aligntype" OR "Pm"
    if "Ii" in kwargs.keys() and ("sigma" in kwargs.keys() or "Pm" in kwargs.keys()):
        Ii = kwargs["Ii"]

        # Validate k
        # if not isinstance(k, (int, float)):
            # raise TypeError(f"k must be a number, got {type(k).__name__}")
        if k < 0:
            raise ValueError(f"k must be non-negative, got {k}")

        # Validate Ii
        # if not isinstance(Ii, (int, float)):
            # raise TypeError(f"'Ii' must be a number, got {type(Ii).__name__}")
        if Ii < 0:
            raise ValueError(f"'Ii' must be non-negative, got {Ii}")

        # Validate sigma/aligntype vs Pm
        if "sigma" in kwargs and "Pm" in kwargs:
            raise ValueError("Provide either 'sigma' and 'aligntype', or 'Pm' — not both.")
        if "sigma" in kwargs and "aligntype" not in kwargs:
            raise ValueError("Option 1 requires both 'sigma' and 'aligntype' when using sigma.")
        if "sigma" in kwargs:
            sigma = kwargs["sigma"]
            aligntype = kwargs["aligntype"]
            if not isinstance(sigma, (int, float)):
                raise TypeError(f"'sigma' must be a number, got {type(sigma).__name__}")
            Pm = calcPm(Ii, sigma, aligntype)
        else:
            Pm = kwargs["Pm"]
            if not isinstance(Pm, dict):
                raise TypeError(f"'Pm' must be a dict mapping substates to probabilities, got {type(Pm).__name__}")

        prefactor = (2*k+1)**(1/2) * (2*Ii + 1)**(1/2)
        total = 0
        Ms = substates(Ii)
        for m in Ms:
            if m not in Pm:
                raise KeyError(f"Substate m={m} expected in 'Pm' but was not found.")
            total += (-1)**(Ii + m) * wigner_3j(Ii, Ii, k, -m, m, 0) * Pm[m]

        return prefactor * total

    # option 2: quantization axis is established by observation of a gamma ray 
    # transition of I0->I1 and a subsequent transition I1->If occurs
    # PARAMETERS: I0, I1, l, lp, delta (mixing ratio)
    elif "I0" in kwargs.keys() and "I1" in kwargs.keys():

        # Check all required keys are present
        required = ["I0", "I1", "l", "lp", "delta"]
        missing = [p for p in required if p not in kwargs]
        if missing:
            raise ValueError(f"Option 2 is missing required parameter(s): {missing}")

        I0 = kwargs["I0"]
        I1 = kwargs["I1"]
        l = kwargs["l"]
        lp = kwargs["lp"]
        delta = kwargs["delta"]

        # Validate types
        for name, val in [("I0", I0), ("I1", I1), ("l", l), ("lp", lp), ("delta", delta)]:
            if not isinstance(val, (int, float)):
                raise TypeError(f"'{name}' must be a number, got {type(val).__name__}")

        # Validate physical constraints
        if I0 < 0:
            raise ValueError(f"'I0' must be non-negative, got {I0}")
        if I1 < 0:
            raise ValueError(f"'I1' must be non-negative, got {I1}")
        if l < 0 or not float(l).is_integer():
            raise ValueError(f"'l' must be a non-negative integer, got {l}")
        if lp < 0 or not float(lp).is_integer():
            raise ValueError(f"'lp' must be a non-negative integer, got {lp}")
        if delta == 0 and l != lp:
            raise ValueError(f"delta=0 implies pure multipolarity, but l={l} != lp={lp}.")
        if 1 + delta**2 == 0:
            raise ZeroDivisionError("Denominator (1 + delta^2) is zero, which is physically invalid.")

        # Note that I0, and I1 are in different positions than typical
        # this is correct and noted properly in https://doi.org/10.1016/S0092-640X(73)80016-6
        return 1.0/(1.0 + delta**2) * (
            get_F(k, l, l, I0, I1)
            + (-1)**(l+lp) * 2 * delta * get_F(k, l, lp, I0, I1)
            + delta**2 * get_F(k, lp, lp, I0, I1)
        )

    else:
        raise ValueError(
            "Invalid or insufficient parameters. Provide either:\n"
            "  Option 1: 'Ii' and ('sigma' + 'aligntype'), or 'Ii' and 'Pm'\n"
            "  Option 2: 'I0', 'I1', 'l', 'lp', 'delta'"
        )