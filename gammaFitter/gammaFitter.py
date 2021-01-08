# gammaFitter.py
# S. Hedges
# 1/8/2021
#
# Description: Code to fit gamma simulation to TOF data for different zero-
# degree backing detector positions.
#
#General approach:
# 1. Load TOF data, applying cut to select higher energy gamma population
#   (onsets are better defined at higher energies). Apply an initial shift this
#   data to this data to bring it into a range to fit.
# 2. Load the simulated data, applying an equivalent energy cut.
# 3. Use emcee to sample different shifts for the simulation data, and offsets
#    for the timing of the simulated depositions, to match the TOF data.Use
#    gaussian (or Gumbel) timing smearing to try to fit any time smearing due to
#    the PMT transit time spread, beam-length, onset smearing, etc.
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
#
# Output:
#   - ROOT file with canvas corresponding to best fit parameters
#   - corner and trace plots in PDF format
import ROOT
import math
import numpy
#EMCEE stuff
import emcee
from multiprocessing import Pool
import matplotlib
from matplotlib import pyplot as plt
import corner
import csv as csvlib
import random
import gc
import os

##############
##User input##
##############
dataFile = "../data/tofData.root"
simFile = "../data/mcnpSims/tofSim_gammas.root"

#Reducing the number of simulation entries speeds of fitting
limitSimToMaxEntries=1 #Flag to decide whether to reduce sim entries
maxSimEntries=2000 #Number of sim entries to reduce to, if above flag set

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

#Guess of fraction of events in gamma peak vs. events in background in fit range
gammaFracGuess=0.85
gammaFracMin=0.75
gammaFracMax=0.95

#Guesses for time shift of simulation to shifted data (in ns)
timeShiftGuess=78.0
timeShiftMin=75.0
timeShiftMax=81.0

#Guesses for time smearing sigma (in ns)
timeSmearGuess=1.2
timeSmearMin=1.1
timeSmearMax=1.3

#EMCEE Parameters
nwalkers=300 #Should be at least 2x ndim!!! I was getting warning messages about
            #initial state not being linearly independent, might have been caused
            #by increasing this too much
#Note this will call the sampler nBurnInSteps*nwalkers times
nBurnInSteps=100
#Note this will call the sampler nSteps*nwalkers times
nSteps=50

multithreaded=1

###################
##RooFit Settings##
###################
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooAbsReal.defaultIntegratorConfig().getConfigSection("RooIntegrator1D").setRealValue("maxSteps",30);
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-8)
ROOT.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-8)

###########################
##Set up RooFit Variables##
###########################
#Range in time to use for PDFs and fits, set so we get a clean gamma population
timeVar=ROOT.RooRealVar("timeVar","timeVar",timeVarMin,timeVarMax)

#Make RooArgSets/Lists once, if created in loop can lead to memory leak
varArgSet=ROOT.RooArgSet(timeVar)
varArgList=ROOT.RooArgList(timeVar)

#############################
##Set up floating variables##
#############################
#Fraction in prompt gamma flash vs. background
gammaFracVar=ROOT.RooRealVar("gammaFracVar","gammaFracVar",gammaFracGuess,gammaFracMin,gammaFracMax)
gammaFracVar.setConstant(1)

#How much time to shift simulation to align with data
timeShiftVar=ROOT.RooRealVar("timeShiftVar","timeShiftVar",timeShiftGuess,timeShiftMin,timeShiftMax)
timeShiftVar.setConstant(1)

#Smearing to apply to simulation to match data
timeSmearVar=ROOT.RooRealVar("timeSmearVar","timeSmearVar",timeSmearGuess,timeSmearMin,timeSmearMax)
timeSmearVar.setConstant(1)

#################
##Load TOF Data##
#################
#Load TOF file for far position, get tree, nEntries
tofDataFile = ROOT.TFile(dataFile,"READ")
zeroDegreeTree = tofDataFile.Get("zeroDegreeTree")
nEntries = zeroDegreeTree.GetEntries()

#Make data set
dataSet=ROOT.RooDataSet("dataSet","dataSet",varArgSet)

#Load data set, shifting data so time always in the correct direction
for entry in range(0,nEntries):
  zeroDegreeTree.GetEntry(entry)
  #Get rid of saturated pulses, pulses outside fit energy
  if zeroDegreeTree.LS_saturated==0 and zeroDegreeTree.LS_integral>ADC_cutOff:
    #If the data is inside the [0,bpmPeriod) BPM range we care about, subtract our hardcoded shift to bring in range
    if zeroDegreeTree.LS_timeToBPM < bpmPeriod and zeroDegreeTree.LS_timeToBPM >= 0:
      time = zeroDegreeTree.LS_timeToBPM
      time -= dataShift
      #If this shifts brings it to negative times, add the BPM period to bring it back in range
      if time < 0:
        time+=bpmPeriod
      #After shifting, only add to PDF if it's the range we're going to fit
      if time >= timeVar.getMin() and time < timeVar.getMax():
        timeVar.setVal(time)
        dataSet.add(varArgSet)
  
print("Loaded TOF data set from file "+dataFile)
  
###############################################
##Load Sim Into Numpy Array For Fast Shifting##
###############################################
gammaSimFile = ROOT.TFile(simFile,"READ")
gammaTree = gammaSimFile.Get("totalEnergyTree")
nEntries = gammaTree.GetEntries()
if nEntries>maxSimEntries and limitSimToMaxEntries==1:
  nEntries=maxSimEntries

#Holds PMT time-smeared data
gammaTimes =[]

#Step through smearing, and adding to list
for entry in range(0,nEntries):
  gammaTree.GetEntry(entry)
  if gammaTree.energy >= keV_cutOff:
    gammaTimes.append(gammaTree.startTime)
  prevHistoryNum=gammaTree.historyNum

#Convert list to numpy array
gammaTimes_arr = numpy.array(gammaTimes)
  
print("Loaded gamma sim from file "+simFile)

################
##Emcee set-up##
################
#Holds the mins and maxes for the emcee call
mins=[]
maxes=[]

labels = ["gammaFrac","timeShift","timeSmearVar"]
mins = [gammaFracVar.getMin(),timeShiftVar.getMin(),timeSmearVar.getMin()]
maxes = [gammaFracVar.getMax(),timeShiftVar.getMax(),timeSmearVar.getMax()]

ndim=len(labels)

#Make into numpy arrays
pos_min = numpy.array(mins)
pos_max = numpy.array(maxes)
#size of each parameter space
psize = pos_max - pos_min
#Generate random values within that space for each walker
pos = [pos_min + psize*numpy.random.rand(ndim) for i in range(nwalkers)]

#Set to 1 after last call to plot best fit restul
plot=0

##############################
########EMCEE FUNCTIONS########
###############################
###Quick shifting of a numpy array###
f = lambda x,shift: x+shift

###Check all values within allowed ranges###
def lnprior(theta):
  #Check if all parameters in range
  allValid=1
  for i in range(0,len(theta)):
    if not mins[i]<theta[i]<maxes[i]:
      allValid=0
    
  if allValid==1:
    return 0.0
  else:
    return -numpy.inf

###Applies time shift to pdfs, forms model using passed in amplitude, and returns pll###
def lnlike(theta):

  #Make a local copy of the amplitude vars so they're not set globally. Set values and make constant
  localGammaFracVar = gammaFracVar.Clone("localGammaFracVar")
  localGammaFracVar.setVal(theta[0])
  localGammaFracVar.setConstant(1)
  
  #Shift numpy array by passed in time shift value
  shiftedTimes=f(gammaTimes_arr,theta[1])
  shiftedTime_arr = numpy.array(shiftedTimes)
 
  #Smear it in time assuming gaussian smearing with passed in sigma
  smearedShiftedTimes=[]
  if smearingType=="gauss":
    for i in range(0,valsToGen):
      smearedShiftedTimes.extend(numpy.random.normal(shiftedTime_arr,theta[2]))
  elif smearingType=="gumbel":
    for i in range(0,valsToGen):
      smearedShiftedTimes.append(numpy.random.gumbel(shiftedTime_arr,theta[2]))
      
  #Make it a flat list
  smearedShiftedTimes_arr=numpy.array(smearedShiftedTimes)
  flatArray=smearedShiftedTimes_arr.flatten()
 
  #Make gamma PDF
  w=numpy.full(flatArray.size,1.)
  h=ROOT.TH1D("h","",nBins,timeVar.getMin(),timeVar.getMax())
  h.FillN(flatArray.size,flatArray,w)
  shiftedSimDataHist=ROOT.RooDataHist("shiftedSimDataHist","shiftedSimDataHist",varArgList,h)
  gammaPdf = ROOT.RooHistPdf("gammaPdf","gammaPdf",varArgSet,shiftedSimDataHist,0)
  h.Delete()
  del h
 
  #Make bgnd PDF
  backgroundPdf = ROOT.RooUniform("backgroundPdf","backgroundPdf",varArgSet)
 
  #Make combined omdel
  pdfList=ROOT.RooArgList(gammaPdf,backgroundPdf)
  ampList=ROOT.RooArgList(localGammaFracVar)
  
  model=ROOT.RooAddPdf("model","model",pdfList,ampList)
  
  ##Generates and returns PLL##
  nll=model.createNLL(dataSet)
  
  #Make it a pll instead of an nll
  pll=-1*nll.getVal()
  
  ##Plotting, only done on last call##
  if plot==1:
    outFile=ROOT.TFile("fitResults_gammas.root","RECREATE")
    c1=ROOT.TCanvas("c1","c1")
    frame = timeVar.frame(timeVar.getMin(),timeVar.getMax(),nBins)
    frame.SetTitle("shift="+str(theta[1])+", smear="+str(theta[2]))
    #Plot source data
    dataSet.plotOn(frame,ROOT.RooFit.Name("Data"),ROOT.RooFit.MarkerColor(1),ROOT.RooFit.FillColor(0))
    model.plotOn(frame,ROOT.RooFit.Name("Model"),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.FillColor(0),ROOT.RooFit.ProjWData(dataSet))
    #Draw
    frame.Draw()
    c1.Modified()
    c1.Update()
    c1.SaveAs("bestFit.png")
    c1.Write()
  #Memory management
  nll.Delete()
  del nll
  model.Delete()
  del model
  backgroundPdf.Delete()
  del backgroundPdf
  gammaPdf.Delete()
  del gammaPdf
  pdfList.Delete()
  del pdfList
  ampList.Delete()
  del ampList
  shiftedSimDataHist.Delete()
  del shiftedSimDataHist
  localGammaFracVar.Delete()
  del localGammaFracVar
  gc.collect

  return pll
  
###Calls our lnlike if passed in values valid###
def lnprob(theta):
  lp = lnprior(theta)
  if not numpy.isfinite(lp):
    return -numpy.inf
  return lp + lnlike(theta)

##############
##Emcee call##
##############

##Multi-threaded##
#Set number of threads to 1 rather than default, will use multiprocessing library for parallelization
if multithreaded==1:
  os.environ["OMP_NUM_THREADS"] = "1"
  with Pool() as pool:

    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool)

    #Burn in only if we're not resuming
    print("Starting burn in...")
    pos, prob, state  = sampler.run_mcmc(pos,nBurnInSteps,progress=True)
    print("Burn-in complete! Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction)))
    sampler.reset()
    pos, prob, state  = sampler.run_mcmc(pos,nSteps,progress=True,store=True)
  
else:
  ##Single threaded##
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
  print("Starting burn in...")
  pos, prob, state  = sampler.run_mcmc(pos, nBurnInSteps, progress=True)
  sampler.reset()
  print("Burn-in complete!")
  pos, prob, state  = sampler.run_mcmc(pos, nSteps, progress=True)

##########################
##Plotting emcee results##
##########################
matplotlib.use('PDF')

#GET THE LL VALUES--FROM GRAYSON'S CODE--Do this first in case plotting fails
samples=sampler.flatchain
lnprobs = sampler.lnprobability[:,:]
flatLnprobs = lnprobs.reshape(-1)
with open("sampler.csv",'a+') as sampleOutFile:
  theWriter = csvlib.writer(sampleOutFile, delimiter=',')
  for sampleLine, llval in zip(samples, flatLnprobs):
    theWriter.writerow(numpy.append(sampleLine,llval))

#MAKE TRACE PLOTS OF EACH PARAMATER
fig = plt.figure(figsize=(10,ndim*2))
gs = fig.add_gridspec(ndim,1)
plt.subplots_adjust(hspace=0.4)
for i in range(0,ndim):
  axes = fig.add_subplot(gs[i,:])
  axes.plot(sampler.chain[:,:,i].T, '-', color='k', alpha=0.3)
  axes.set_title(labels[i])
plt.savefig("traceplots.pdf")

#CORNER PLOT HERE
#Make extents for corner plots
extents=[]
for i in range(0,len(mins)):
  entry=[]
  entry.append(mins[i])
  entry.append(maxes[i])
  extents.append(entry)
samples=sampler.flatchain
fig = corner.corner(samples, labels=labels, range=extents, quantiles=[0.16,0.5,0.84],show_titles=True,title_kwargs={'fontsize':12})
fig.savefig("corner.pdf")

#CALCULATE QUANTILES HERE
bestFitValues=[]
for i in range(ndim):
  mcmc = numpy.percentile(sampler.chain[:,:,i],[16, 50, 84])
  q = numpy.diff(mcmc)
  print(labels[i]+": "+str(mcmc[1])+"+"+str(q[0])+" -"+str(q[1])+"\n")
  bestFitValues.append(mcmc[1])

#Plot with best values
plot=1
lnlike(bestFitValues)

#Print to screen--do this at end because sometimes fails!
print("Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction))+"\n")
print("Mean autocorrelation time: {0:.3f} steps".format(numpy.mean(sampler.get_autocorr_time(c=1,quiet=True))))
