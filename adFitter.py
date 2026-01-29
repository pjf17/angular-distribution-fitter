import sys
import warnings
import numpy as np
from ROOT import TF1, TLine, TLatex, TGraphErrors, TLegend, TCanvas, gStyle

sys.path.append("modules")
from measurement import Measurement
from angular_analyze import AngularDistribution, adRatios, readSourceFile, readDataFile, readStatsFile, doFit

NUC = {'12Ex': 0.335}

def psuedoH0(Xdata,Yerr):
    #generate pseudo data according to the null hypothesis (isotropic distribution)
    Ydata = np.random.normal(1.0,Yerr) 
    return Ydata

def psuedoH1(Xdata,Yerr,AngDist):
    #generate pseudo data according to the fit hypothesis
    Ydata = np.random.normal(AngDist.evalFitW(Xdata),Yerr) 
    return Ydata

if __name__ == "__main__":
    if len(sys.argv) < 9:
        print("Not enough arguments, please enter the following arguments:")
        print("<source file> <data file> <Beta> <Algn Type (P/O)> <E> <Ji> <Jf> <Mix> <(optional) dual option>")
        print("The last argument can either be:") 
        print("> number of pseudo experiments to run to obtain parameter distributions")
        print("OR")
        print("> .stats file after analyzing pseudodata to add the p-values to the plot legend")
        sys.exit(1)

    #read in the data
    theta, srce = readSourceFile(sys.argv[1])
    data = readDataFile(sys.argv[2])
    peak = int(sys.argv[5])

    # read in the optional last argument
    nExp = 0 #number of pseudo data experiments
    H0_pval = H1_pval = None #get p-values for dchi2 and chi2
    if len(sys.argv) == 10:
        if sys.argv[9].isnumeric(): nExp = int(sys.argv[9])
        else: H0_pval, H1_pval = readStatsFile(sys.argv[9])

    #get the beta, check for nucleus code
    beta = 0
    plttitle = ""
    if sys.argv[3] in NUC:
        beta = NUC[sys.argv[3]]
        plttitle = "^{" + sys.argv[3][:2] + "}" + sys.argv[3][2:]  
        plttitle = plttitle + ", E_{#gamma} = " + sys.argv[5] 
        plttitle = plttitle + ", #beta = " + f"{beta:.3f}"
    else:
        beta = float(sys.argv[3])
        plttitle = "#beta = " + f"{beta:.3f}"

    #initialize angular distribution properties
    AD = AngularDistribution(sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[4],beta,theta)

    #get the angular distribution ratios
    ratioData = adRatios(peak,beta,theta,data[peak],srce)
    AD.calcFlatChi2(ratioData)
    rdV  = []
    rdE = []
    for rd in ratioData:
        rdV.append(rd.value)
        rdE.append(rd.error)
    
    theta = np.array(theta)
    rdV = np.array(rdV)
    rdE = np.array(rdE)

    #use root functions to do the fit and plot the result
    #the fit returns the same result as angular_analyze.doFit
    #but I prefer to do the actual plotting using root functions
    if nExp == 0:
        thetaError = np.zeros(len(theta)) #TGraphErrors needs a x error array as well, set all to zero

        #make the scatter plot
        ADpoints = TGraphErrors(len(theta),theta,rdV,thetaError,rdE)
        ADpoints.SetMarkerStyle(8)
        ADpoints.SetMarkerSize(1.5)
        ADpoints.SetLineWidth(2)

        #do the fit
        wrapper = lambda x, par: AD.TF1W(x,par)
        f1 = TF1("fleg", wrapper, min(theta), max(theta), 2)
        f1.SetLineWidth(3)
        f1.SetParameter(0,0.6)
        f1.SetParLimits(0,0.0,1000.0)
        if AD.fixMix:
            f1.FixParameter(1,AD.delta.value)
        else:
            f1.SetParameter(1,0.0)
            f1.SetParLimits(1,-10.0,10.0)
        fitres = ADpoints.Fit(f1,"EX0 N0R S")
        AD.getFitResult(fitres.Chi2(),fitres.Ndf(),Measurement(fitres.Parameter(0),fitres.ParError(0)),Measurement(fitres.Parameter(1),fitres.ParError(1)))
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
        if H1_pval is not None: legend.AddEntry(0,"p(#chi^{2})= "+f"{H1_pval:.2f}","")
        legend.AddEntry(flat,"Iso","l")
        if H0_pval is not None: legend.AddEntry(0,"p(#Delta#chi^{2})= "+f"{H0_pval:.2f}","")
        
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

    # Don't plot the fit, but instead run psuedo experiments to get the uncertainties on the
    # fit parameters
    if nExp > 0:
        # angular_analyze.doFit is faster than the root fitter, so we use it do iterate over nExp
        print("Real data fit results:")
        doFit(AD,theta,rdV,rdE)
        AD.Print()
        print(f"Running {nExp} pseudo data simulations")

        #format the Ji, Jf strings for better filename format
        fileJi = 0
        fileJf = 0
        if isinstance(AD.Ii,int):
            fileJi = AD.Ii
            fileJf = AD.If
        else:
            fileJi = f"{AD.Ii}".replace('/','-')
            fileJf = f"{AD.If}".replace('/','-')
        
        outFile = f"pseudo_data_{peak}_{fileJi}_{fileJf}.dat"
        with open(outFile, "w") as file:
            file.write("h0_dchi2,h0_chi2,h0_fltchi2,h1_dchi2,h1_chi2,h1_fltchi2,align,pol,sigma,delta")
            for k in AD.a.keys():
                file.write(f",a{k}")
            file.write("\n")
            file.write(f"{AD.fltChi2-AD.chi2},{AD.chi2},{AD.fltChi2},{AD.fltChi2-AD.chi2},{AD.chi2},{AD.fltChi2},{AD.calcAlign()},{AD.calcPol()},{AD.sigma.value},{AD.delta.value}")
            for k in AD.a.keys():
                file.write(f",{AD.a[k]}")
            file.write("\n")

            #initialize the distributions under the two hypotheses
            AD_H0 = AngularDistribution(sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[4],beta,theta)
            AD_H1 = AngularDistribution(sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[4],beta,theta)

            #run the loop
            for i in range(nExp):
                #H0
                psV = psuedoH0(theta,rdE)
                AD_H0.calcFlatChi2(psV,rdE)
                doFit(AD_H0,theta,psV,rdE)

                #H1
                psV = psuedoH1(theta,rdE,AD)
                AD_H1.calcFlatChi2(psV,rdE)
                doFit(AD_H1,theta,psV,rdE)

                file.write(f"{AD_H0.fltChi2-AD_H0.chi2},{AD_H0.chi2},{AD_H0.fltChi2},{AD_H1.fltChi2-AD_H1.chi2},{AD_H1.chi2},{AD_H1.fltChi2},{AD_H1.calcAlign()},{AD_H1.calcPol()},{AD_H1.sigma.value},{AD_H1.delta.value}")
                for k in AD_H1.a.keys():
                    file.write(f",{AD_H1.a[k]}")
                file.write("\n")

                AD_H0.resetStats()
                AD_H1.resetStats()
                print(f"\r{i+1}/{nExp}",end='')
            
            print(f"\nFinished, data saved to {outFile}")
