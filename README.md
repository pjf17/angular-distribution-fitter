# Angular Distribution Fitter
This python script will plot and fit angular distribution data to obtain the alignment and degree of polarization percentages. It uses [SymPy](https://www.sympy.org/en/index.html), [Numpy](https://numpy.org), [SciPy](https://scipy.org), and [PyROOT](https://root.cern/manual/python/) (tested using with ROOT version 6.28.04). 

![](output.svg)
## Running
The figure above was generated using the following command:
`python3 -i adFitter.py source.dat example_data.dat 12Ex Obl 2311 4 2 0.0`

To generate pseudo data run the same command with the number of experiments to run:
`python3 -i adFitter.py source.dat example_data.dat 12Ex Obl 2311 4 2 0.0 500`

To make histograms for all fit parameter and derived parameter distributions, use:
`python3 makeDistHists.py <your psuedo_data.dat file>`

`makeDistHists.py` also generates a `.stats` file that contains hypothesis test and best fit p-values that can be added to the angular distribution fit legend. This can be plotted using:
`python3 -i adFitter.py source.dat example_data.dat 12Ex Obl 2311 4 2 0.0 <your matching .stats file>`

`makeDistHists.py` also generates a `.root` file containing all histograms of the parameter distributions that can be inspected in the root browser

## Input Files
The data file is formatted where each line contains all of the angular data for a transition in a single line. The format is: `<peak> <theta bin 1 counts> <theta bin 1 error> ...` where any number of theta bins > 1 can be used. See `example_data.dat` 

The source file contains the 152Eu source information. The first column is the center of each theta bin. The number of theta bins for the source file must match the number for the data file. The second column is the number of counts at the corresponding theta for the 1408 peak. The third column is the error on the counts. See `source.dat`

## Formalism
The angular distribution $W(\theta)$ is the distribution of gamma-ray intensity, which requires knowing $\epsilon(E,\theta_i)$, which is the efficiency of observing the lab frame gamma-ray energy in the polar angle cut at $\theta_i$. This efficiency is calculated according to the following. The total GRETINA singles efficiency for a gamma ray that is Doppler shifted according to angle $\theta_i$ is $\epsilon\left(E \gamma (1-\beta\cos\theta_i)\right)$ where $\epsilon(E)$ is given by equation \ref{eq:efficiency-curve}. The total efficiency at the Doppler-shifted energy is then scaled according to the fractional angular coverage in GRETINA of the angle bin. The fractional angular coverage can be found by sorting source data according to the same polar angle cuts and obtaining the full energy peak counts of a transition within each angle cut. The 1408-keV peak from ${}^{152}\mathrm{Eu}$ was used to obtain the fractional angular coverage for GRETINA. The efficiency $\epsilon(E,\theta)$ can be expressed as:

$\epsilon(E,\theta_i) = \epsilon(E\gamma(1-\beta\cos\theta_i)) \frac{N^{1408}_i}{\sum_j{N^{1408}_j}}$

where $N^{1408}_i$ is the number of 1408-keV ${}^{152}\mathrm{Eu}$ source counts within the angle bin $\theta_i$. This is implemented in the function `adRatios()`.