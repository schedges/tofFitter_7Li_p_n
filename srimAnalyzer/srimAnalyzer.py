# srimAnalyzer.py
# S. Hedges
# 1/8/2021
#
# Description: Uses simulations of protons incident on a LiF film at various
# depths to produce 0-degree neutron energy distributions. To convert from
# proton to neutron energy, interpolates values from Liskien and Paulsen
# (https://doi.org/10.1016/0092-640X(75)90004-2) correcting for changes to the
# neutron production cross section with energy.
#
# Usage: Specify SRIM files, input target thickness, and step size. Output is a
# ROOT file with proton and neutron energy distributions, and MCNP-like sources
# printed to the screen.
#
# Notes:
#   - Near the 7Li(p,n) production threshold, Lee & Zhou
#     (https://doi.org/10.1016/S0168-583X(99)00026-9) data should be used
#     instead.
#   - The code does not correct for any change in the proton flux throughout
#     the target. This is a small effect for thin targets, but may be important
#     for thicker or gas targets.
#   - There could be issues with the MCNP output if there are gaps in the
#     neutron energy distribution. Please check for gaps in the PoliMi output
#     lines greater than the supplied bin size, if using.
#
import ROOT
import math

##############
##User input##
##############
#Specify target thickness
targetThickness = 750 #nm

#Specify SRIM files, incident proton energies
energies = [2680,2690,2700]
filenames = ["../data/srimSims/"+str(i)+"keV_25nm_475nm.txt" for i in energies]

#Output bin size in keV
binSize=2 #keV

#Output filename
outputFile=ROOT.TFile("srimOutputs.root","RECREATE")

#Ignore any bins with less than 1% of neutrons to simplify source definition
threshold_fraction=0.01 #fraction

################################
##Liskien & Paulsen input data##
################################
#Proton energies (keV) from Liskine & Paulsen
e_p = [1950,2000,2050,2150,2200,2250,2300,2350,2400,2450,2500,2600,2700,2800,
       2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,
       4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,
       5700,5800,5900,6000,6100,6200,6300,6400,6500,6600,6700,6800,6900]
#0-degree neutron energies corresponding to supplied proton energies (keV)
e_n = [165,230,291,350,407,463,518,573,627,680,733,891,996,1099,
       1203,1304,1408,1511,1613,1715,1816,1918,2019,2121,2222,2323,2424,2525,
       2626,2727,2828,2929,3030,3130,3231,3332,3432,3533,3633,3734,3835,3935,
       4035,4136,4236,4336,4437,4538,4638,4738,4838,4939,5039,5139,5239]
#0-degree neutron x-sections corresponding to supplied proton energies (mb/sr)
xs_n = [58.8,37.8,27.2,44.6,88.5,145.,149.,124.,104.,89.3,78.7,65.6,57.2,53.5,
        50.9,48.8,47.1,45.7,44.4,43.2,42.1,41.1,40.3,39.6,38.8,38.3,38.7,40.1,
        42.6,45.6,48.9,52.8,27.8,62.8,67.3,39.9,67.8,63.5,58.9,55.5,50.7,46.9,
        43.7,40.0,36.8,33.6,31.1,28.8,26.8,24.9,23.1,21.3,19.9,18.7,17.6]
  
###############################
##Make TGraphs to interpolate##
###############################
#Make graphs of liskien-paulsen proton energies vs. n energies, vs. x-sections
liskienNeutronEnergyGraph=ROOT.TGraph()
liskienNeutronXSectionGraph=ROOT.TGraph()
nPoints=0

#All we care about is relative x-sections, so scale based on first x-section
xs_n = [i/xs_n[0] for i in xs_n] #Make relative x-sections

#Fill graphs
for i in range(0,len(e_p)):
  liskienNeutronEnergyGraph.SetPoint(nPoints,e_p[i],e_n[i])
  liskienNeutronXSectionGraph.SetPoint(nPoints,e_p[i],xs_n[i])
  nPoints+=1

######################
##Load up SRIM files##
######################
dataStartLine=13 #Skip SRIM header info from transit file

#Step through files and generate an array of proton energies
protonEnergyLists=[]
for i,filename in enumerate(filenames):
  protonEnergyLists.append([])

for i,filename in enumerate(filenames):
  #Open srim file
  f = open(filename,"r")
  
  #Keep track of line number
  lineNo=1
  
  #Step through lines
  for line in f:
    #Skip header
    if lineNo>=dataStartLine and not "===" in line:
      #Strip newline character
      line=line.strip("\n")
      #Split based on what space
      lineParts=line.split()
      #Get proton energy in keV, it starts in eV
      protonEnergy=float(lineParts[3])/1000.
      #Add to histogram
      protonEnergyLists[i].append(protonEnergy)
    lineNo+=1
  f.close()

################################################
##Convert from proton energy to neutron energy##
################################################
neutronHist_min=0
neutronHist_max=6000
neutronHist_bins=int((neutronHist_max-neutronHist_min)/binSize)

#Now step through and convert from proton energy to neutron energy
neutronEnergyHists=[]
for i,protonEnergyList in enumerate(protonEnergyLists):
  #Make neutron hist
  neutronEnergyHists.append(ROOT.TH1D("neutronEnergyHist_"+str(energies[i])+"keV",";Neutron energy (keV);Counts",neutronHist_bins,neutronHist_min,neutronHist_max))
  
  #Step through proton energies
  for protonEnergy in protonEnergyList:
  
      #Get neutron energy, weighting factor
      neutron_energy=liskienNeutronEnergyGraph.Eval(protonEnergy)
      weightingFactor = liskienNeutronXSectionGraph.Eval(protonEnergy)
      
      #Get neutron energy bin
      neutron_energy_bin=neutronEnergyHists[i].FindBin(neutron_energy)
      
      #Add the proton bin contents to any content already in neutron energy bin
      neutronEnergyHists[i].SetBinContent(neutron_energy_bin,neutronEnergyHists[i].GetBinContent(neutron_energy_bin)+weightingFactor)

########################################
##Scale and write neutron energy hists##
########################################
for i,neutronEnergyHist in enumerate(neutronEnergyHists):
  print(str(energies[i])+":")

  #Get integral
  integral=neutronEnergyHist.Integral()
  #Normalize so integral is 1.
  neutronEnergyHist.Scale(1./integral)
  
  #Generate MCNP-like output lines. Start with bin edges of histogram
  siLine = "si2 H "
    
  for ibin in range(1,neutronEnergyHist.GetNbinsX()+1):
    if neutronEnergyHist.GetBinContent(ibin)>=threshold_fraction:
    
      #Get bin center in keV
      binCenter=neutronEnergyHist.GetBinCenter(ibin)
      
      #Convert to MeV
      binCenter_MeV = binCenter/1000
      
      #Get left edge of bin for PoliMi source
      #Get the bin size in MeV
      binSize_MeV = float(binSize)/1000.
      #Subtract 1/2 bin size from bin center to get bin edge
      binEdge = binCenter_MeV - binSize_MeV/2.
      #Write line
      siLine += "{0:.3f} ".format(binEdge)
        
  #Polimi requires an upper bin edge
  lineParts=siLine.split()
  lastEnergy=float(lineParts[-1])
  finalBinEdge=lastEnergy+binSize_MeV
  siLine+=str(finalBinEdge)
  
  #Print to screen
  print(siLine)
  
  #Now supply bin probabilities, Everything to right of first bin assumed 0
  spLine = "sp2 0 "
  for ibin in range(0,neutronEnergyHist.GetNbinsX()+1):
    if neutronEnergyHist.GetBinContent(ibin)>=threshold_fraction:
      spLine += "{0:.3f} ".format(neutronEnergyHist.GetBinContent(ibin))
     
  #Print to screen
  print(spLine+"\n")
  
  #Save histgrams
  outputFile.cd()
  neutronEnergyHist.Write()
