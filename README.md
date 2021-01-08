# tofFitter_7Li_p_n
Code to fit TOF distributions from 7Li(p,n) reactions. Code has three components:

# srimAnalyzer.py
Takes SRIM outputs of proton energy distributions different distances into a thin LiF target, and combines to get a realistic energy distribution of protons in the whole target. Converts to 0-degree neutron energy distributions using energy and cross-section data from Liskien & Paulsen (https://doi.org/10.1016/0092-640X(75)90004-2). Outputs histograms, as well as MCNP-compatible histogram source.

# gammaFitter.py
Fits simulations of gammas originating on the LiF target to TOF data to determine a shift to apply to simulations, as well as a time-resolution smearing function. 

# neutronFitter.py
Using the best-fit shift and smearing from the gamma fitter, and simulations of neutrons generated on the LiF target with the distributions pulled from srimAnalyzer.py, creates a RooMomentMorph as a function of incident proton energy. Fits to find the best fit incident energy, and then generates the corresponding neutron distribution. 
