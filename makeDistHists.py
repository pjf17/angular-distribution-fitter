import sys
import os
import math
import pandas as pd
import numpy as np
from collections import namedtuple
from ROOT import TFile, TCanvas, TH1F, TLine

def getBounds(hist,blo=15.9,bhi=84.1):
    # for handling non-symmetric error bars
    # get the 16th and 84th percentiles, which is the same as a
    # 1 sigma confidence region as a gaussian (68% confidence interval)
    total = hist.Integral()
    nbins = hist.GetNbinsX()

    boundLo = boundHi = 0
    for b in range(2,nbins+1):
        frac = hist.Integral(1,b)*1.0/total
        if frac > blo/100:
            boundLo = hist.GetBinLowEdge(b)
            break

    for b in range(nbins-1,1,-1):
        frac = hist.Integral(1,b)*1.0/total
        if frac < bhi/100:
            boundHi = hist.GetBinLowEdge(b+1)
            break

    return boundLo, boundHi

def setLineProperties(l,Width=None,Style=None,Color=None):
    if Width is not None:
        l.SetLineWidth(Width)
    if Style is not None:
        l.SetLineStyle(Style)
    if Color is not None:
        l.SetLineColor(Color)

def printLatex(name,center,down,up):
    # output central values and bounds already formatted in LaTeX 
    # for easy copy-pasting into overleaf
    hiBound = loBound = None
    if center < 0.0:
        loBound = abs(center-up)
        hiBound = abs(center-down)
    else:
        hiBound = abs(center-up)
        loBound = abs(center-down)
    output = "$" + f"{center:.3f}" + "^{+" + f"{hiBound:.3f}" +  "}_{-" + f"{loBound:.3f}" + "}" + "$"
    print(name, output)

#settings for each histogram for the parameters
HIST_INIT = {"h0_dchi2": (30000,-5,25), "h1_chi2": (50000,0,50), "align": (1000,0,100), "pol": (1000,0,100), 
                "sigma": (100,0,2.0), "delta": (200,-10,10), "a2": (4000,-0.5,0.5), "a4": (4000,-0.5,0.5)}

#choose which histograms to plot, all get saved to a .root file for future inspection
DISPLAY_DATA = ["h0_dchi2","h1_chi2","align","pol"]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_histogram.py <filename.dat>")
        sys.exit(1)

    filename = sys.argv[1]
    
    #if asymmetry + error is given, also create & draw the pol sensitivity (Q) hist
    asym = None
    asymerr = None
    if len(sys.argv) == 4:
        asym = float(sys.argv[2])
        asymerr = float(sys.argv[3])
        DISPLAY_DATA.append("Q")
        HIST_INIT.update({"Q": (200,0,2.0)})

    try:
        # Load full data
        df_full = pd.read_csv(filename)

        # Extract first row as best-fit values
        best_fit_row = df_full.iloc[0]
        BestFit = namedtuple("BestFit", df_full.columns)
        best_fit = BestFit(*best_fit_row)

        # Remaining data (exclude first row)
        df = df_full.iloc[1:]

        #set the output files
        outname = os.path.splitext(filename)[0] + ".root" #save the histograms here
        statsname = os.path.splitext(filename)[0] + ".stats" #for adding stats to AD plot
        fout = TFile(outname,"RECREATE")

        #prep the canvas to draw 
        c = TCanvas("Canvas", "Canvas", 800, 600)
        Ndisp = len(DISPLAY_DATA)
        ny = 2
        nx = (Ndisp + ny - 1) // ny
        c.Divide(nx,ny)

        hists = []
        iDisp = 0
        with open(statsname,"w") as file:
            for column_name in HIST_INIT.keys():
                doDraw = column_name in DISPLAY_DATA 

                # Create histogram
                hists.append( TH1F(f"{column_name}", f"{column_name}", HIST_INIT[column_name][0], HIST_INIT[column_name][1], HIST_INIT[column_name][2]) )

                # Fill the histogram from column values if parameter is not Q
                if column_name is not "Q":
                    values = (pd.to_numeric(df[column_name], errors='coerce').dropna()).tolist()
                    for val in values:
                        hists[-1].Fill(val)

                # Fill the Q histogram, where Q = asym/pol and asym follows a gaus distribution
                else:
                    polvalues = (pd.to_numeric(df["pol"], errors='coerce').dropna()).tolist()
                    for val in polvalues:
                        hists[-1].Fill(np.random.normal(asym,asymerr)/val*100)

                # Set best fit line
                bfx = getattr(best_fit, column_name) if column_name is not "Q" else asym/(getattr(best_fit, "pol")/100)
                line = TLine(bfx,0,bfx,hists[-1].GetMaximum())
                hists[-1].GetListOfFunctions().Add(line)
                setLineProperties(line,2,1,2)

                #get upper and lower bounds if not a chi2
                if "chi2" not in column_name:
                    blo, bhi = getBounds(hists[-1])
                    llo = TLine(blo,0,blo,hists[-1].GetMaximum())
                    lhi = TLine(bhi,0,bhi,hists[-1].GetMaximum())
                    setLineProperties(llo,2,2,1)
                    setLineProperties(lhi,2,2,1)
                    hists[-1].GetListOfFunctions().Add(llo)
                    hists[-1].GetListOfFunctions().Add(lhi)
                    printLatex(column_name,bfx,blo,bhi)
                    file.write(f"{column_name}:{bfx} {blo} {bhi}\n")
                
                #get the p values on the hypothesis test
                if (column_name is "h0_dchi2") or (column_name is "h1_chi2"):
                    # h0_dchi2 hypothesis test (chi2h0 - chi2h1) assuming my data followed h0
                    # what is the probability that the difference in chi2 between h1 and h0 is larger
                    # or equal to the value of my experiment, if Pval is high, then the null
                    # hypothesis stands (sad), if Pval is low then we reject null hypothesis (yay)

                    # h1_chi2, assuming the data follows h1, what percent of pseudo_data events are
                    # worse than or equal to my value, if Pval is high, I have a good fit
                    bestFitValBin = hists[-1].FindBin(getattr(best_fit, column_name))
                    print(f"Pval {column_name}:",hists[-1].Integral(bestFitValBin,hists[-1].GetNbinsX())*1.0/hists[-1].Integral())
                    file.write(f"{column_name}:{hists[-1].Integral(bestFitValBin,hists[-1].GetNbinsX())*1.0/hists[-1].Integral()}\n")

                #save hist to file
                hists[-1].Write()

                # Draw
                if doDraw:
                    c.cd(iDisp+1)
                    hists[-1].Draw()
                    iDisp += 1

    except Exception as e:
        print("An error occurred:", e)
