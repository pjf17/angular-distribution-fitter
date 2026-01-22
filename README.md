# Angular Distribution Fitter
This python script will plot and fit angular distribution data to obtain the alignment and degree of polarization percentages. It uses [SymPy](https://www.sympy.org/en/index.html), [Numpy](https://numpy.org), [SciPy](https://scipy.org), and [PyROOT](https://root.cern/manual/python/) (tested using with ROOT version 6.28.04). 

![](output.svg)

The figure above was generated using the following command:
`python3 -i fitdata.py source.dat example_data.dat 12Ex Obl 2311 4 2 0.0`