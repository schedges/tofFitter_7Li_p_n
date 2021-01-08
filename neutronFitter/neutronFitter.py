# neutronFitter.py
# S. Hedges
# 1/8/2021
#
# Description: Using supplied shift and smearing parameters from the
# gammaFitter.py code, and simulations of neutrons energy distributions
# (obtained via SRIM and Liskien & Paulsen calculations) originating on target,
# fits the incident proton energy. Based on that best fit value, creates the
# best fit neutron energy distribution.
#
# Makes some assumptions on the format of the simulation and data files, based
# on the included files.
#   Simulation tree name: totalEnergyTree
#   Simulation energy branch name: energy (keV)
#   Simulation event start time name: startTime (ns)
#
#   Data tree name: zeroDegreeTree
#   Data tree time to BPM branch: LS_timeToBPM (ns)
#   Data tree is saturated branch: LS_isSaturated
#   Data tree energy branch: LS_integral (ADC)
#   Data tree PSD branch: LS_psd
#
# Notes:
#   - Settings should be copied over from the gamma fitting code.
import ROOT
import random
import os
import numpy

###########################
##Sim and data file names##
###########################
dataFile = "../data/tofData.root"
#Sorted sim energies and sim files
energies = [2680,2690,2700]
simFiles = ["../data/mcnpSims/tofSim_2680keV_n.root","../data/mcnpSims/tofSim_2690keV_n.root","../data/mcnpSims/tofSim_2700keV_n.root"]

# Range to select a clean neutron population from data, after dataShift applied.
# Also want a relatively flat backgrounds over this range.
fitRangeMin = 120 #ns
fitRangeMax = 155 #ns

#Lower psd cut-off for neutrons
psd_cutOff_lo = 0.24
#Higher psd cut-off for neutrons
psd_cutOff_hi =0.55

#File with srim distributions--from srimAnalyzer.py"
srimDistributionsFile=ROOT.TFile("../srimAnalyzer/srimOutputs.root","READ")

#Reducing the number of simulation entries speeds of fitting
limitSimToMaxEntries=0 #Flag to decide whether to reduce sim entries
maxSimEntries=10000 #Number of sim entries to reduce to, if above flag set

############################################
##Input best fit shift from gammaFitter.py##
############################################
simShift = 78.878 #ns
#Input best fit smearing sigma from gammaFitter.py
smearingSigma = 1.246 #ns

###########################################
##Copy these over from gamma fitting code##
###########################################
#Number of smeared values to generate for each simulated event
valsToGen=20

#Which type of smearing function to use
smearingType = "gauss" #"gauss" or "gumbel"

#Apply an initial shift to the data so gammas and neutrons don't wrap around the
#BPM period.
dataShift = 150 #ns
  
#Nominal time between BPM pulses
bpmPeriod=400 #ns
  
#ADC->keV conversion: E = ADC/c + Eo
# Used to apply an energy cut in the simulation to select similar neutron
# population (energy resolution not taken into account)
c = 18.79 #ADC/keV
Eo = 5 #keV
ADC_cutOff = 2000 #ADC
keV_cutOff = ADC_cutOff/c + Eo

#Gamma fitting region after dataShift has been applied
timeVarMin=0 #ns
timeVarMax=120 #ns

#Binning--how many bins/ns
binsPerNS = 4
nBins = int(binsPerNS*(timeVarMax-timeVarMin))

#################
##RooFit set-up##
#################
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL) #DANGER

#Make RooRealVar
timeVar=ROOT.RooRealVar("timeVar","timeVar",fitRangeMin,fitRangeMax)
timeBinning = ROOT.RooBinning(nBins,timeVar.getMin(),timeVar.getMax(),"timeBinning")
timeVar.setBinning(timeBinning)

#Make arg list and arg set
varArgSet=ROOT.RooArgSet(timeVar)
varArgList=ROOT.RooArgList(timeVar)

#Load in TOF data set
dataSet=ROOT.RooDataSet("dataSet","dataSet",varArgSet)
tofFile = ROOT.TFile("tofData.root","READ")
tree=tofFile.Get("zeroDegreeTree")
nEntries = tree.GetEntries()
for entry in range(0,nEntries):
  tree.GetEntry(entry)
  if tree.LS_saturated==0 and tree.LS_integral >= ADC_cutOff and tree.LS_psd >= psd_cutOff_lo and tree.LS_psd<psd_cutOff_hi:
    time = tree.LS_timeToBPM-dataShift
    if time < 0:
      time+=bpmPeriod
    if time >= timeVar.getMin() and time < timeVar.getMax():
      timeVar.setVal(time)
      dataSet.add(varArgSet)
      
print("Loaded TOF data")

#For holding generated
hists=[]
dataHists=[]
histPdfs=[]
pdfs = ROOT.RooArgList()
paramVec = ROOT.TVectorD(len(energies))

#Generate pdfs:
for j,simFile in enumerate(simFiles):
  print(energies[j],simFile)
  inpFile = ROOT.TFile(simFile,"READ")
  polimiTree = inpFile.Get("totalEnergyTree")
  nEntries = polimiTree.GetEntries()
  if nEntries>maxSimEntries and limitSimToMaxEntries==1:
    nEntries=maxSimEntries
    
  #Make a TH1
  hists.append(ROOT.TH1D(str(energies[j])+"keV",str(energies[j])+"keV",nBins,timeVar.getMin(),timeVar.getMax()))
  
  #For checking if we need to add it to the histogram
  prevHistoryNum=0
  
  for entry in range(0,nEntries):
    polimiTree.GetEntry(entry)
    
    #If it's a new history and the projectile is a neutron and energy > cut off
    if polimiTree.quenchedEnergy >= keV_cutOff:
    
      #Shift the timestamp
      timestamp = polimiTree.startTime + simShift
      
      #Apply smearing
      for i in range(0,valsToGen):
        if smearingType=="gauss":
          valToFill =  numpy.random.normal(timestamp,smearingSigma)
        else:
          valToFill =  numpy.random.gumbel(timestamp,smearingSigma)
        
        #After smearing, if we're outside of the range bring it back in.
        if valToFill > bpmPeriod:
          valToFill-=bpmPeriod
        elif valToFill < 0:
          valToFill += bpmPeriod
          
        #Fill our histogram
        if valToFill >= timeVar.getMin() and valToFill < timeVar.getMax():
          hists[-1].Fill(valToFill)

  #Convert hist to RooDataHist -> RooHistPdf
  dataHists.append(ROOT.RooDataHist(hists[-1].GetName()+"_DataHist",hists[-1].GetName()+"_DataHist",varArgList,hists[-1]))
  histPdfs.append(ROOT.RooHistPdf(hists[-1].GetName()+"_HistPdf",hists[-1].GetName()+"_HistPdf",varArgSet,dataHists[-1],1))
  
  #Add to rooarglist
  pdfs.add(histPdfs[-1])
  paramVec[j] = energies[j]

  #Close file
  inpFile.Close()
       
print("Loaded sims")

#Make morphing parameter
protonEnergyVar = ROOT.RooRealVar("protonEnergyVar","protonEnergyVar",energies[0],energies[0],energies[-1])
#What kind of morph
setting = ROOT.RooMomentMorph.Linear
#The morph pdf
morphPdf = ROOT.RooMomentMorph("morphPdf","morphPdf",protonEnergyVar,varArgList,pdfs,paramVec,setting)
neutronFracVar=ROOT.RooRealVar("neutronFracVar","neutronFracVar",0.5,0,1.0)

#The background PDF
backgroundPdf = ROOT.RooUniform("backgroundPdf","backgroundPdf",varArgSet)

#The combined PDF
model=ROOT.RooAddPdf("model","model",ROOT.RooArgList(morphPdf,backgroundPdf),ROOT.RooArgList(neutronFracVar))

##Generates and returns PLL##
res=model.fitTo(dataSet,
  ROOT.RooFit.Extended(1),
  ROOT.RooFit.Save(True),
  ROOT.RooFit.NumCPU(4)
)
  
c1=ROOT.TCanvas("c1","c1")

frame = timeVar.frame(timeVar.getMin(),timeVar.getMax(),nBins)
bestEnergy = protonEnergyVar.getVal()
err = protonEnergyVar.getError()
frame.SetTitle(str(bestEnergy)+"+-"+str(err)+" protons")

#Plot source data
dataSet.plotOn(frame,ROOT.RooFit.Name("Data"),ROOT.RooFit.MarkerColor(1),ROOT.RooFit.FillColor(0))

dataHist=ROOT.RooDataHist("dataHist","dataHist",varArgSet,dataSet)
tofDataHist = dataHist.createHistogram("tofDataHist",timeVar,ROOT.RooFit.Binning(nBins))
integral=tofDataHist.Integral()
tofDataHist.Scale(1./integral)
tofSimHist = model.createHistogram("tofSimHist",timeVar,ROOT.RooFit.Binning(nBins))

model.plotOn(frame,ROOT.RooFit.Name("Model"),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.FillColor(0))
model.plotOn(frame,ROOT.RooFit.Name("Neutrons"),ROOT.RooFit.Components("morphPdf"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(dataSet))
model.plotOn(frame,ROOT.RooFit.Name("Background"),ROOT.RooFit.Components("backgroundPdf"),ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGray),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(dataSet))
  
#Draw
c1.cd()
frame.Draw()

leg = ROOT.TLegend(0.65,0.65,0.95,0.92);
dataObj = frame.findObject("Data");
modelObj = frame.findObject("Model");
bgndObj = frame.findObject("Background");
neutronObj = frame.findObject("Neutrons");
leg.AddEntry(dataObj,"Data","P")
leg.AddEntry(modelObj,"Model","L")
leg.AddEntry(bgndObj,"Background","L")
leg.AddEntry(neutronObj,"Neutrons","L")

leg.Draw("same")
c1.Modified()
c1.Update()

energyHists=[]
energyDataHists=[]
energyHistPdfs=[]
energyPdfs = ROOT.RooArgList()
  
neutronEnergyVar=ROOT.RooRealVar("neutronEnergyVar","neutronEnergyVar",0,6000)
neutronEnergyVarBins=3000
energyBinning = ROOT.RooBinning(neutronEnergyVarBins,neutronEnergyVar.getMin(),neutronEnergyVar.getMax(),"energyBinning")
neutronEnergyVar.setBinning(energyBinning)
neutronEnergyVarArgList=ROOT.RooArgList(neutronEnergyVar)
neutronEnergyVarArgSet=ROOT.RooArgSet(neutronEnergyVar)

#Make histograms
for i in range(0,len(energies)):
  name = "neutronEnergyHist_"+str(energies[i])+"keV"
  energyHists.append(srimDistributionsFile.Get(name))
  energyDataHists.append(ROOT.RooDataHist(str(energies[i])+"_EnergyDataHist",str(energies[i])+"_EnergyDataHist",neutronEnergyVarArgList,energyHists[i]))
  energyHistPdfs.append(ROOT.RooHistPdf(str(energies[i])+"_EnergyHistPdf",str(energies[i])+"_EnergyHistPdf",neutronEnergyVarArgSet,energyDataHists[i],0)) #No interpolation
  energyPdfs.add(energyHistPdfs[i])
  
#The energy morph pdf
energyMorphPdf = ROOT.RooMomentMorph("energyMorphPdf","energyMorphPdf",protonEnergyVar,neutronEnergyVarArgList,energyPdfs,paramVec,setting)
protonEnergyVar.setVal(bestEnergy)
protonEnergyVar.setError(err)

outFile=ROOT.TFile("fitResults_neutrons.root","RECREATE")
c1.SetName("tofDataAndFit")
c1.Write()
energyHist = energyMorphPdf.createHistogram("energyHist",neutronEnergyVar,ROOT.RooFit.Binning(neutronEnergyVarBins))
energyHist.SetName("neutronEnergyDistribution_"+str(bestEnergy)+"keV")
energyHist.Write()
tofDataHist.SetName("tofData_neutrons")
tofDataHist.Write()
tofSimHist.SetName("tofSim_neutrons")
tofSimHist.Write()
outFile.Close()
