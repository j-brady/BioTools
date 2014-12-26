BioTools
========
Some scripts for analysis of protein sequences

	def example(test)
		return test
## Dependencies
It is necessary to install the Numpy, MatPlotLib and BioPython for some of these scripts to work.

	pip install numpy matplotlib BioPython
## Installation

Download/clone the package. In the rootdir of the package run the setup.py script.

	python setup.py install


## data_tables.py
Contains hydropathy tables.

## WebTools
Scripts for incorporation into a Django web application

NMR Tools
=========
Some scripts for dealing with NMR related data
## nus.py
Scripts for converting Bruker NUS sampling schedules to Omega compatible format and vice-versa.

	from NMR_Tools import nus
	
	nus.Omega2Bruk('nuslist1','nuslist2') 

will output a bruker format nus list for use with hmsIST

	nus.Bruk2Omega('nuslist')

will convert bruker nuslist to Omega format


## Tract
Scripts for analysis of <sup>15</sup>N-TRACT data. Fitting and effective correlation time (&tau;<sub>c</sub>) along with estimation of apparent molecular weight.  

###Example usage:

	from NMR_Tools.tract import tract_tauc, Mw_from_Tauc
        
	tau_c = tract_tauc.tauc(h_field,two_eta_xy)
        
h_field = field strength in Hz, two_eta_xy = cross-correlated relaxation rate (2&eta;<sub>xy</sub>) (s<sup>-1</sup>) calculated from <sup>15</sup>N-TRACT data. 	
	
	mw    = Mw_from_Tauc(r_w,eta,T,tau_c,rho)

r_w = hydration radius (should be between 1.6 - 3.2 &#8491;), eta = viscosity, tau_c= &tau;<sub>c</sub> (s), rho = 0.73e-6 m<sup>3</sup>g<sup>-1</sup>.
