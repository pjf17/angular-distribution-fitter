import numpy as np
from scipy.special import legendre
from scipy.special import lpmv
from scipy.optimize import curve_fit
from measurement import Measurement
import angular_formalism as adf

#================================================================================
#functions for reading in the data, source, and stats files, these each contain different
#types of info, so they have slightly different file formats
#================================================================================
def readDataFile(filename):
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

def readSourceFile(filename):
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
    # read in values from the stats files so they can be added to the TLegend
    try:
        data_dict = {}
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split(':')
                data_dict[parts[0]] = parts[1]

            pH0 = float(data_dict["h0_dchi2"])
            pH1 = float(data_dict["h1_chi2"])
            # align = (float(data_dict["align"].split()[1]),float(data_dict["align"].split()[2]))
            # pol = (float(data_dict["pol"].split()[1]),float(data_dict["pol"].split()[2]))

        return pH0, pH1 #, align, pol
    
    except FileNotFoundError:
        raise FileNotFoundError("File '{:f}' not found".format(filename))

#================================================================================
#For calculating the angular distribution ratios
#================================================================================
def calcEff(E, xtalNum=48):
    #get GRETINA efficiency based on Dirk's NIM calibration curve for however many xtals you have
    return 4.532*xtalNum/32*pow(E+100,-0.621)

def dopShift(E, beta, theta):
    gamma = 1/np.sqrt(1-beta*beta)
    return E/(gamma*(1-beta*np.cos(theta*np.pi/180)))

def adRatios(energy, beta, theta, data, src):
    #get the angular distribution ratios, normalized to an "isotropic" distribution 
    #in the lab frame 
    ratio = []

    sourceSum = Measurement(0,0)
    dataSum = Measurement(0,0)
    N = int(len(data))
    for i in range(N):
        sourceSum = sourceSum + src[i]
        dataSum = dataSum + data[i]
        temp = data[i] / (calcEff(dopShift(energy,beta,theta[i]),43))
        ratio.append(temp / src[i].value)

    ratioSum = Measurement(0,0)
    for i in range(N):
        ratio[i] = ratio[i] * sourceSum.value
        ratio[i] = ratio[i] / ( (1-beta*beta)/(beta*np.cos(theta[i]*np.pi/180)-1)**2 )
        ratioSum = ratioSum + ratio[i]
    
    for i in range(N):
        ratio[i] = ratio[i] / ratioSum.value * N
        effEff = calcEff(dopShift(energy,beta,theta[i]),43)*src[i].value/sourceSum.value *( (1-beta*beta)/(beta*np.cos(theta[i]*np.pi/180)-1)**2 )

    return ratio

#================================================================================
# non-ROOT fitting function
#================================================================================
def doFit(AngDist,theta,ratioVal,ratioErr):
    my_fit = None
    my_bounds = None

    if AngDist.fixMix:
        my_fit = lambda x, s: AngDist.W(x,[s,AngDist.delta.value])
        my_bounds = ([0],[np.inf])
    else:
        my_fit = lambda x, s, d: AngDist.W(x,[s,d])
        my_bounds = ([0,-10.0],[np.inf,10.0])
    
    #do the fit
    popt, pcov = curve_fit(my_fit, theta, ratioVal, sigma=ratioErr, bounds=my_bounds, maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    popt = np.array(popt)

    #calc the chi2
    ndf = len(theta) - len(popt)
    chi2 = 0
    for i in range(len(theta)):
        chi2 += ((ratioVal[i]-my_fit(theta[i], *popt))/ratioErr[i] )**2

    sigma = Measurement(popt[0],perr[0])
    if AngDist.fixMix:
        AngDist.getFitResult(chi2,ndf,sigma)
    else:
        delta = Measurement(popt[1],perr[1])
        AngDist.getFitResult(chi2,ndf,sigma,delta)

#================================================================================
class AngularDistribution:
    # class to handle Angular Distribution fitting
    def __init__(self, Ii, If, delta, alignment, beta, theta):
        self.Ii = self.readSpin(Ii) #initial spin state
        self.If = self.readSpin(If) #final spin state
        self.beta = beta #beam velocity

        self.align = None #alignment type
        if alignment == "Pro": # True for prolate, False for oblate
            self.align = True 
        else:
            self.align = False
        
        self.delta = None #mixing ratio
        self.fixMix = False #fix the mixing ratio
        if delta is not '-':
            self.fixMix = True
            self.delta = Measurement(float(delta),0.0)

        #get the lowest possible multipolarity of the transition
        self.L = abs(self.Ii - self.If)
        if self.L is 0: self.L = 1 #gamma transitions must be at least L=1

        #cache the legendre calculations for the reference points
        tmpTh = np.array(theta)
        tmpTh = (np.cos(tmpTh/180*np.pi) - self.beta)/(1 - self.beta*np.cos(tmpTh/180*np.pi))
        self.legendreThetaPoints = {k: legendre(k)(tmpTh) for k in range(0,2*(self.L+2),2)}
        
        #variables that will be used later
        self.sigma = None #alignment width parameter
        self.A = None #A coefficients
        self.Ap = None #A' coefficients 
        self.B = None #B coeffcieients
    
    def readSpin(self,I):
        #from user input, determine whether spin is half integer or integer
        if '/' in I:
            return S(I)
        else:
            return int(I)
    
    def __str__(self):
        output = f"{self.chi2},{self.ndf},{self.fltChi2},{self.fltNdf},{self.sigma.value},{self.delta.value},{self.calcAlign()},{self.calcPol()}"
        for k in self.a.keys():
            output += f",{self.a[k]}"
        return output
    
    def resetStats(self):
        self.chi2 = None
        self.ndf = None
        self.A = None
        self.Ap = None
        self.B = None
        self.a = None
        self.fltChi2 = None
        self.fltNdf = None

    def getFitResult(self, chi2, ndf, sigma, delta=None):
        self.chi2 = chi2
        self.ndf = ndf
        self.sigma = sigma
        if delta is not None:
            self.delta = delta

        self.A = {x: float(adf.A_k(x,self.delta.value,self.L,self.L+1,self.If,self.Ii)) for x in range(0,2*(self.L+2),2)}
        self.Ap = {x: float(adf.Ap_k(x,self.delta.value,self.L,self.L+1,self.If,self.Ii)) for x in range(0,2*(self.L+2),2)}
        self.B = {x: float(adf.B_k(x,self.Ii,self.sigma.value,self.align)) for x in range(0,2*(self.L+2),2)}
        self.a = {x: self.A[x]*self.B[x] for x in range(0,2*(self.L+2),2)}
    
    def calcPol(self,lo = 47,hi = 107,n=100):
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
        return float(self.B[2]/adf.B_k(2,self.Ii,0.0,self.align))*100

    def calcFlatChi2(self,ratio,err=None):
        #get chi2 under assumption that the distribution is constant =1
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

    def TF1W(self, x, par):
        #the angular distribution fit function for ROOT Tf1
        #the fit parameters enter into the coefficients, par[0] = sigma, par[1] = delta
        A = {k: float(adf.A_k(k,par[1],self.L,self.L+1,self.If,self.Ii))*float(adf.B_k(k,self.Ii,par[0],self.align)) for k in range(0,2*(self.L+2),2)}
        cosCm = (np.cos(x[0]*np.pi/180) - self.beta)/(1 - self.beta*np.cos(x[0]*np.pi/180))
        
        norm = 0
        w = 0
        for k in A.keys():
            w += A[k]*legendre(k)(cosCm)
            for tp in self.legendreThetaPoints[k]:
                norm += A[k]*tp
        
        return w/norm * len(self.legendreThetaPoints[0])
    
    def W(self, x, par):
        #the angular distribution fit function 
        #the fit parameters enter into the coefficients, par[0] = sigma, par[1] = delta
        A = {k: float(adf.A_k(k,par[1],self.L,self.L+1,self.If,self.Ii))*float(adf.B_k(k,self.Ii,par[0],self.align)) for k in range(0,2*(self.L+2),2)}
        cosCm = (np.cos(x*np.pi/180) - self.beta)/(1 - self.beta*np.cos(x*np.pi/180))
        
        norm = 0
        w = 0
        for k in A.keys():
            w += A[k]*legendre(k)(cosCm)
            for tp in self.legendreThetaPoints[k]:
                norm += A[k]*tp
        
        return w/norm * len(self.legendreThetaPoints[0])
    
    def evalFitW(self, x):
        #evaluate the best fit function
        cosCm = (np.cos(x*np.pi/180) - self.beta)/(1 - self.beta*np.cos(x*np.pi/180))

        norm = 0
        w = 0
        for k in self.a.keys():
            w += self.a[k]*legendre(k)(cosCm)
            for tp in self.legendreThetaPoints[k]:
                norm += self.a[k]*tp
        
        return w/norm * len(self.legendreThetaPoints[0])