from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from sympy.physics.quantum.cg import CG
from sympy import S
import sympy as sp
import math
import numpy as np
from scipy.special import legendre
from scipy.special import lpmv
import sys
import warnings
from ROOT import TF1, TLine, TLatex, TGraphErrors, TLegend, TCanvas, gStyle

class Measurement:
    #class for handling error propagation
    def __init__(self, value, error):
        self.value = float(value)
        self.error = float(abs(error))  # Ensure error is positive
    
    def __str__(self):
        return f"{self.value:.3} ± {self.error:.3}"
    
    def __repr__(self):
        return f"Measurement({self.value}, {self.error})"
    
    def __add__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value + other.value
        new_error = math.sqrt(self.error**2 + other.error**2)
        return Measurement(new_value, new_error)
    
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value - other.value
        new_error = math.sqrt(self.error**2 + other.error**2)
        return Measurement(new_value, new_error)
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        new_value = self.value * other.value
        rel_error1 = self.error / abs(self.value) if self.value != 0 else 0
        rel_error2 = other.error / abs(other.value) if other.value != 0 else 0
        new_error = abs(new_value) * math.sqrt(rel_error1**2 + rel_error2**2)
        return Measurement(new_value, new_error)
    
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        if other.value == 0:
            raise ValueError("Division by zero")
        new_value = self.value / other.value
        rel_error1 = self.error / abs(self.value) if self.value != 0 else 0
        rel_error2 = other.error / abs(other.value) if other.value != 0 else 0
        new_error = abs(new_value) * math.sqrt(rel_error1**2 + rel_error2**2)
        return Measurement(new_value, new_error)
    
    # To allow operations with numbers on the left (e.g., 2 * Measurement)
    def __radd__(self, other):
        return self.__add__(other)
    
    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        return other.__sub__(self)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            other = Measurement(other, 0)
        return other.__truediv__(self)

#================================================================================
#multipolarity coefficients/functions
#================================================================================
def get_kap(k, l, lp):
    if (k > l + lp or k == 0):
        return 0
    prefactor = -1*np.sqrt(np.math.factorial(k-2)/np.math.factorial(k+2))
    cgnum = CG(l,1,lp,1,k,2)
    cgden = CG(l,1,lp,-1,k,0)
    return prefactor*float(cgnum.doit())/float(cgden.doit())

def get_F(k, l, lp, If, Ii):
    prefactor = (-1)**(If+Ii+1) * (2*k+1)**(1/2) * (2*l+1)**(1/2) * (2*lp+1)**(1/2) * (2*Ii+1)**(1/2)
    return prefactor * wigner_3j(l,lp,k,1,-1,0) * wigner_6j(l,lp,k,Ii,Ii,If)

def get_u(k, l, Ii, If):
    prefactor = (-1)**(Ii+If+l+k) * ((2*Ii+1)*(2*If+1))**(1/2)
    return prefactor * wigner_6j(Ii,Ii,k,If,If,l)

def get_U(k, delta, l, lp, Ii, If):
    return (get_u(k,l,Ii,If) + delta**2*get_u(k,lp,Ii,If))/(1+delta**2)

def get_A(k, delta, l, lp, If, Ii):
    return 1.0/(1.0 + delta**2) * (get_F(k,l,l,If,Ii) + 2*delta*get_F(k,l,lp,If,Ii) +  delta**2*get_F(k,lp,lp,If,Ii))

def get_Ap(k, delta, l, lp, If, Ii):
    return 1.0/(1.0 + delta**2) * (-1.0*get_kap(k,l,l)*get_F(k,l,l,If,Ii) + 2*delta*get_kap(k,l,lp)*get_F(k,l,lp,If,Ii) +  delta**2*get_kap(k,lp,lp)*get_F(k,lp,lp,If,Ii))

#================================================================================
#alignment functions
#================================================================================
def getSubstates(I):
    Ms = []
    mm = -I
    for i in range(int(2*I)+1):
        Ms.append(mm)
        mm += 1
    
    return Ms

def getPm(I,sigma,algntype):
    # calculates the substate population distribution under differen alignment type assumptions
    Pm = {}
    norm = 0
    subs = getSubstates(I)
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

def get_B(k, Ii, sigma, algntype,Pm=None):
    #the alignment coefficient
    prefactor = (2*k+1)**(1/2) * (2*Ii + 1)**(1/2)
    sum = 0
    Ms = getSubstates(Ii)

    if (Pm == None):
        Pm = getPm(Ii,sigma,algntype)

    for m in Ms:
        sum += (-1)**(Ii + m)*wigner_3j(Ii,Ii,k,-m,m,0)*Pm[m]
    
    return prefactor*sum

#================================================================================
#functions for reading in the data, source, and stats files
#these each contain different types of info, so they have slightly different file formats
#================================================================================
def readDatFile(filename):
    data = {}
    
    try:
        with open(filename, 'r') as file:
            for line_number, line in enumerate(file, 1):
                # Skip empty lines
                if line.strip() == "":
                    continue
                    
                # Split the line into parts and remove extra whitespace
                parts = line.split()

                peak = int(parts[0])
                nDataPoints = int((len(parts)-1)/2)
                points = []
                for i in range(nDataPoints):
                    val = float(parts[i*2+1])
                    err = float(parts[i*2+2])
                    points.append(Measurement(val,err));

                data.update({peak: points})
                    
        return data
    
    except FileNotFoundError:
        raise FileNotFoundError("File '{:s}' not found".format(filename))

def readSrcFile(filename):
    data = []
    theta = []
    
    try:
        with open(filename, 'r') as file:
            for line_number, line in enumerate(file, 1):
                # Skip empty lines
                if line.strip() == "":
                    continue
                    
                # Split the line into parts and remove extra whitespace
                parts = line.split()

                theta.append(float(parts[0]))
                data.append(Measurement(float(parts[1]),float(parts[2])))
                
                    
        return theta, data
    
    except FileNotFoundError:
        raise FileNotFoundError("File '{:f}' not found".format(filename))

def readStatsFile(filename):
    try:
        data_dict = {}
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split(':')
                data_dict[parts[0]] = parts[1]

            pH0 = float(data_dict["h0_dchi2"])
            pH1 = float(data_dict["h1_chi2"])
            align = (float(data_dict["align"].split()[1]),float(data_dict["align"].split()[2]))
            pol = (float(data_dict["pol"].split()[1]),float(data_dict["pol"].split()[2]))

        return pH0, pH1, align, pol
    
    except FileNotFoundError:
        raise FileNotFoundError("File '{:f}' not found".format(filename))

def calcEff(E):
    #get GRETINA efficiency
    return 4.532*43.0/32*pow(E+100,-0.621)

def dopShift(E, beta, theta):
    gamma = 1/np.sqrt(1-beta*beta)
    return E/(gamma*(1-beta*np.cos(theta*np.pi/180)))

def readSpin(I):
    #from user input, determine whether spin is half integer or integer
    if '/' in I:
        return S(I)
    else:
        return int(I)

def getRatios(energy, beta, theta, data, src):
    #get the angular distribution ratios, normalized to an "isotropic" distribution 
    #in the lab frame 
    ratio = []

    sourceSum = Measurement(0,0)
    dataSum = Measurement(0,0)
    N = int(len(data))
    for i in range(N):
        sourceSum = sourceSum + src[i]
        dataSum = dataSum + data[i]
        temp = data[i] / (calcEff(dopShift(energy,beta,theta[i])))
        ratio.append(temp / src[i].value)

    ratioSum = Measurement(0,0)
    for i in range(N):
        ratio[i] = ratio[i] * sourceSum.value
        ratio[i] = ratio[i] / ( (1-beta*beta)/(beta*np.cos(theta[i]*np.pi/180)-1)**2 )
        ratioSum = ratioSum + ratio[i]
    
    for i in range(N):
        ratio[i] = ratio[i] / ratioSum.value * N
        effEff = calcEff(dopShift(energy,beta,theta[i]))*src[i].value/sourceSum.value *( (1-beta*beta)/(beta*np.cos(theta[i]*np.pi/180)-1)**2 )

    return ratio

class AngularDistribution:
    # class to handle Angular Distribution fitting
    def __init__(self, Ii, If, delta, alignment, beta, theta):
        self.Ii = Ii #initial spin state
        self.If = If #final spin state
        self.beta = beta #beam velocity

        self.align = None #alignment type
        if alignment == "Pro": # True for prolate, False for oblate
            self.align = True 
        else:
            self.align = False
        
        self.fixMix = False #fix the mixing ratio
        if delta is not '-':
            self.fixMix = True
            self.delta = float(delta)

        #get the lowest possible multipolarity of the transition
        self.L = abs(self.Ii - self.If)
        if self.L is 0: self.L = 1 #gamma transitions must be at least L=1

        #cache the legendre calculations for the reference points
        tmpTh = np.array(theta)
        tmpTh = (np.cos(tmpTh/180*np.pi) - self.beta)/(1 - self.beta*np.cos(tmpTh/180*np.pi))
        self.legendreThetaPoints = {k: legendre(k)(tmpTh) for k in range(0,2*(self.L+2),2)}
    
    def __str__(self):
        output = f"{self.chi2},{self.ndf},{self.fltChi2},{self.fltNdf},{self.sigma.value},{self.delta.value},{self.calcAlign()},{self.calcPol()}"
        for k in self.a.keys():
            output += f",{self.a[k]}"
        return output

    def getFitResult(self, fitres):
        self.chi2 = fitres.Chi2()
        self.ndf = fitres.Ndf()
        self.sigma = Measurement(fitres.Parameter(0),fitres.ParError(0))
        self.delta = Measurement(fitres.Parameter(1),fitres.ParError(1))
        self.A = {x: float(get_A(x,self.delta.value,self.L,self.L+1,self.If,self.Ii)) for x in range(0,2*(self.L+2),2)}
        self.Ap = {x: float(get_Ap(x,self.delta.value,self.L,self.L+1,self.If,self.Ii)) for x in range(0,2*(self.L+2),2)}
        self.B = {x: float(get_B(x,self.Ii,self.sigma.value,self.align)) for x in range(0,2*(self.L+2),2)}
        self.a = {x: self.A[x]*self.B[x] for x in range(0,2*(self.L+2),2)}
    
    def calcPol(self,lo = 48,hi = 115,n=100, printTerms=False):
        #calculate the degree of Polarization percentage
        theta = np.linspace(lo*np.pi/180,hi*np.pi/180,n)
        legP = {}
        assocP = {}
        legSum = 0
        assocSum = 0
        N = len(theta)
        for k in self.A.keys():
            tmpLeg = 0
            tmpAssoc = 0
            for t in theta:
                cosTheta = (np.cos(t) - self.beta)/(1 - self.beta*np.cos(t))
                tmpLeg += legendre(k)(cosTheta)
                tmpAssoc += lpmv(2,k,cosTheta) 

            legP[k] = tmpLeg/N
            assocP[k] = tmpAssoc/N

            legSum += self.A[k]*self.B[k]*legP[k]
            assocSum += self.Ap[k]*self.B[k]*assocP[k]

        total = abs(assocSum/legSum)

        return total*100

    def calcAlign(self):
        #return the Alignment percentage
        return float(self.B[2]/get_B(2,self.Ii,0.0,self.align))*100

    def calcFlatChi2(self,ratio,err=None):
        self.fltChi2 = 0
        self.fltNdf = len(ratio)-1
        if err is None:
            for r in ratio:
                self.fltChi2 += ((r.value - 1.0)/(r.error))**2
        else: 
            for i in range(len(ratio)):
                self.fltChi2 += ((ratio[i] - 1.0)/(err[i]))**2

    def Print(self):
        print("=============== FIT RESULTS ===============")
        print(f"chi2/ndf (flat): {self.fltChi2:.2f}/{self.fltNdf } = {self.fltChi2/self.fltNdf :.2f}")
        print(f"chi2/ndf  (fit): {self.chi2:.2f}/{self.ndf} = {self.chi2/self.ndf:.2f}")
        print("Sig/J:",self.sigma)
        print("Delta:",self.delta)
        print("=============== CALCULATIONS ==============")
        print("Ji -> Jf\tAlgn\tPol",end="")
        for k in self.a.keys():
            print(f"\ta{k}",end="")
        print("")
        print(f" {self.Ii} -> {self.If} \t{self.calcAlign():3.1f}%\t{self.calcPol():3.1f}%",end="")
        for k in self.a.keys():
            print(f"\t{self.a[k]:.3f}",end="")
        print("")

    def W(self, x, par):
        #the angular distribution fit function
        #the fit parameters enter into the coefficients, par[0] = sigma, par[1] = delta
        A = {k: float(get_A(k,par[1],self.L,self.L+1,self.If,self.Ii))*float(get_B(k,self.Ii,par[0],self.align)) for k in range(0,2*(self.L+2),2)}
        
        cosCm = (np.cos(x[0]*np.pi/180) - self.beta)/(1 - self.beta*np.cos(x[0]*np.pi/180))
        norm = 0
        w = 0
        for k in A.keys():
            w += A[k]*legendre(k)(cosCm)
            for tp in self.legendreThetaPoints[k]:
                norm += A[k]*tp
        
        return w/norm * len(self.legendreThetaPoints[0])

NUC = {'12Ex': 0.335}

if __name__ == "__main__":
    if len(sys.argv) < 9:
        print("<source file> <data file> <Beta> <Algn Type (P/O)> <E> <Ji> <Jf> <Mix>")
    else: 
        #read in the data files
        theta, srce = readSrcFile(sys.argv[1])
        data = readDatFile(sys.argv[2])
        peak = int(sys.argv[5])

        #get the beta, check for nucleus code, set up plot title
        plttitle = ""
        beta = 0
        if sys.argv[3] in NUC:
            beta = NUC[sys.argv[3]]
            plttitle = "^{" + sys.argv[3][:2] + "}" + sys.argv[3][2:]  
            plttitle = plttitle + ", #beta = " + f"{beta:.3f}"
            plttitle = plttitle + ", E_{#gamma} = " + sys.argv[5] 
        else:
            beta = float(sys.argv[3])

        #initialize angular distribution properties
        AD = AngularDistribution(readSpin(sys.argv[6]),readSpin(sys.argv[7]),sys.argv[8],sys.argv[4],beta,theta)

        #prep the data
        ratioData = getRatios(peak,beta,theta,data[peak],srce)
        AD.calcFlatChi2(ratioData)
        rdV  = []
        rdE = []
        for rd in ratioData:
            rdV.append(rd.value)
            rdE.append(rd.error)

        rdV = np.array(rdV)
        rdE = np.array(rdE)
        thx = np.array(theta)
        thxE = np.zeros(len(thx))
        ADpoints = TGraphErrors(len(thx),thx,rdV,thxE,rdE)
        ADpoints.SetMarkerStyle(8)
        ADpoints.SetMarkerSize(1.5)
        ADpoints.SetLineWidth(2)

        #do the fit
        wrapper = lambda x, par: AD.W(x,par)
        f1 = TF1("fleg", wrapper, min(thx), max(thx), 2)
        f1.SetLineWidth(3)
        f1.SetParameter(0,0.6)
        f1.SetParLimits(0,0.0,1000.0)
        if AD.fixMix:
            f1.FixParameter(1,AD.delta)
        else:
            f1.SetParameter(1,0.0)
            f1.SetParLimits(1,-10.0,10.0)
        AD.getFitResult(ADpoints.Fit(f1,"EX0 N0R S"))
        AD.Print()

        #set up plot
        ymax = 1.5
        ymin = 0.7
        gStyle.SetPadLeftMargin(0.12)
        gStyle.SetPadRightMargin(0.04)
        gStyle.SetPadTopMargin(0.02)
        legend = TLegend(0.65,0.78,0.96,0.98)
        canv = TCanvas("canv","canv",600,600)

        #draw the frame
        haxis = canv.DrawFrame(0,ymin,180,ymax)

        #configure axes
        haxis.GetXaxis().SetTitle("#theta [#circ]")
        haxis.GetXaxis().SetTitleSize(0.055)
        haxis.GetXaxis().SetTitleOffset(0.9)
        haxis.GetXaxis().SetLabelSize(0.042)
        haxis.GetYaxis().SetTitle("W(#theta)/W(#theta)_{iso}")
        haxis.GetYaxis().SetTitleSize(0.05)
        haxis.GetYaxis().SetTitleOffset(1.08)
        haxis.GetYaxis().SetLabelSize(0.042)
        
        #configure the flat function
        flat = TLine(0,1.0,180,1.0)
        flat.SetLineColor(1)
        flat.SetLineStyle(2)
        flat.SetLineWidth(2)

        #configure the fit function to be drawn
        fitFunc = TF1("fleg", wrapper, 0,180, 2)
        fitFunc.SetParameters(f1.GetParameters())
        legend.AddEntry(fitFunc,f"{AD.Ii} #rightarrow {AD.If}","l")
        legend.AddEntry(flat,"Iso","l")
        
        #draw Plot title
        latex = TLatex()
        latex.SetTextSize(0.048)
        latex.SetTextFont(42)
        latex.DrawLatex(5,(ymax-ymin)*0.92 + ymin,plttitle)

        #draw everything
        ADpoints.Draw("P same")
        fitFunc.Draw("same")
        legend.Draw("same")
        flat.Draw("same")
        