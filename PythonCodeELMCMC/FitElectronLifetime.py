################################
### Will use emcee for fitting
################################
import ElectronLifetimeTrend
from ElectronLifetimeTrend import *

#import MyHistorianLib
#from MyHistorianLib import GetUnixTimeFromTimeStamp

import FormPars

import Tools
from Tools import *

import ROOT
from ROOT import TCanvas 
from ROOT import TFile
from ROOT import TTree
from ROOT import TChain
from ROOT import TCut
from ROOT import TH1
from ROOT import TH1D
from ROOT import TH2
from ROOT import TH2D
from ROOT import TGraph
from ROOT import TGraphErrors
from ROOT import TF1
from ROOT import TLatex
from ROOT import TLine
from ROOT import TPad

import numpy as np
from scipy.optimize import minimize
from numpy.linalg import inv

import time
import copy
import click
import pickle
import emcee


StartingTimeFit = time.time()
if len(sys.argv)<3:
    print("======== Syntax ==========")
    print("python FitElectronLifetime.py .......")
    print("<SC file>")
    print("<Output pickle result>")
    print("<Fit data txt>")
    print("<Rn Fit data txt>")
    print("<number of walkers>")
    print("<number of iteractions>")
    print("<(optional)Input pickle pre-result>")
    exit()

HistorianFile = sys.argv[1]
FitOutput = sys.argv[2]
S1ExponentialConstant = 2040.6 # us. 

print('\nFitting Electron Lifetime between ' + FormPars.GetMinTimeStamp() + ' and ' + FormPars.GetMaxTimeStamp() + '\n')

# setting the parameters
MinUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMinTimeStamp())
MaxUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMaxTimeStamp())
default_pars = FormPars.GetDefaultPars()

# initial parameters
x0, x0_steps = FormPars.GetInitialParametersMCMC()

# The main Light yield Trend
pElectronLifetimeTrend = MyElectronLifetimeTrend(HistorianFile, MinUnixTime, MaxUnixTime, default_pars)

# Need to load the S2/S1 elife data
ElectronLifetimeDataFile = sys.argv[3]
f1 = open(ElectronLifetimeDataFile)
lines1 = f1.readlines()
f1.close()

# Need to also load Rn data
RnElectronLifetimeDataFile = sys.argv[4]
f2 = open(RnElectronLifetimeDataFile)
lines2 = f2.readlines()
f2.close()


# pre-walking
nwalkers = int(sys.argv[5])
niterations = int(sys.argv[6])
PreWalkingPickleFilename = "NoneExist"
if len(sys.argv)>7:
    PreWalkingPickleFilename = sys.argv[7]

ElectronLifetimeData = {}
gEL = TGraphErrors()
#############################
# fill in the data
#############################
# electron lifetime
UnixTimes = []
Values = []
ValueErrors = []
for i, line in enumerate(lines1):
    contents = line[:-1].split("\t\t")
    #print(contents[0])
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    # correct the e-life by the S1 term
    value = value*S1ExponentialConstant / (S1ExponentialConstant - value)
    value_err = value*np.sqrt( np.power(value_err/value, 2.0)+np.power(value_err/ (S1ExponentialConstant - value), 2.) )
    if value_err<=0:
        continue
    UnixTimes.append(unixtime)
    Values.append(value)
    ValueErrors.append(value_err)
    gEL.SetPoint(i, unixtime, value)
    gEL.SetPointError(i, unixtime_err, value_err)


## Load the Rn data
RnUnixtimes = []
RnUnixtimeErrors = []
RnELifeValues = []
RnELifeValueErrors = []

for line2 in lines2:
    if line2[0] == '#':
        continue
    contents = line2[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    #if was calculated using pax_v6.2.0 or younger
#    if unixtime < 1478000000 or (unixtime > 1484900000 and unixtime < 1486100000):
    if unixtime < 1478000000:
        value, value_err = CorrectForPaxVersion(value, value_err)
    RnUnixtimes.append(unixtime)
    RnUnixtimeErrors.append(unixtime_err)
    RnELifeValues.append(value)
    RnELifeValueErrors.append(value_err)

RnMinUnixtime = np.min(RnUnixtimes)
# Remove values that has unixtime larger than the first one in RnUnixTimes
# And then extend the list with Rn lists
CutOffID = 0
for i, unixtime in enumerate(UnixTimes):
    if unixtime > RnMinUnixtime:
        CutOffID = i
        break
UnixTimes = UnixTimes[0:CutOffID]
Values = Values[0:CutOffID]
ValueErrors = ValueErrors[0:CutOffID]
UnixTimes.extend(RnUnixtimes)
Values.extend(RnELifeValues)
ValueErrors.extend(RnELifeValueErrors)



ElectronLifetimeData['UnixTimes'] = UnixTimes
ElectronLifetimeData['Values'] = Values
ElectronLifetimeData['ValueErrors'] = ValueErrors
LastPointUnixTime = UnixTimes[len(UnixTimes)-1]

for unixtime in ElectronLifetimeData['UnixTimes']:
    print( unixtime )

# log likelihood function for mcmc fitting
def LnLikeData():
    UnixTimes = ElectronLifetimeData['UnixTimes']
    Values = ElectronLifetimeData['Values']
    ValueErrors = ElectronLifetimeData['ValueErrors']
    LnL = 0.
    for unixtime, value, value_err in zip(UnixTimes, Values, ValueErrors):
        expected = pElectronLifetimeTrend.GetElectronLifetime(unixtime)
        LnL += -0.5* np.power((value - expected)/value_err, 2.)
    return LnL

def LnLike(x):
    start_time = time.time()
    pars, IfOutOfBoundary = FormPars.FormPars(x)
    if IfOutOfBoundary:
#        print("Out of boundary!")
        return -np.inf
    global pElectronLifetimeTrend
    pElectronLifetimeTrend.SetParameters(pars)
#    print("==== The parameters: =====")
#    print(pElectronLifetimeTrend.GetParameters())
    LnL = LnLikeData()
#    print(" LnL = "+str(LnL))
#    print("=====================")
#    print("Takes: "+str(time.time()-start_time) + " sec")
    return LnL

####################
# start MCMC fitting
####################
ndim = len(x0)
# randomize the walkers
p0 = []
IfPickleSuccessful = False
if not PreWalkingPickleFilename=="NoneExist":
    PreWalkingData = pickle.load(open(PreWalkingPickleFilename, 'rb'))
    nwalkers = PreWalkingData['nwalkers']
    chain = PreWalkingData['chain']
    acceptance_fractions = PreWalkingData['acceptance_fraction']
    PreviousSample = chain[0][-1]
    for Sample, acceptance_fraction in zip(ReshapeChain(chain, ndim, nwalkers, PreWalkingData['niterations'])[-1], acceptance_fractions):
        if acceptance_fraction<0.10:
            p0.append(PreviousSample)
            continue
        p0.append(Sample)
        PreviousSample = Sample
    IfPickleSuccessful = True
elif not IfPickleSuccessful:
    for i in range(nwalkers):
#        p0.append([np.random.normal(a,b) for a, b in zip(x0, x0_steps)])
        p0.append([np.random.uniform(a-2*b,a+2*b) for a, b in zip(x0, x0_steps)])

#sampler = emcee.EnsembleSampler(nwalkers, ndim, LnLike, threads=28)
# from Matt's github https://github.com/mdanthony17/emcee
sampler = emcee.DESampler(nwalkers, ndim, LnLike, threads=28)
#sampler = emcee.DESampler(nwalkers, ndim, LnLike)

#for results in sampler.sample(p0, iterations=niterations):
#    pass
#with sampler.sample(p0=p0, iterations=niterations) as mcmc_sampler:

with click.progressbar(sampler.sample(p0=p0, iterations=niterations), length=niterations) as mcmc_sampler:
    for i,results in enumerate(mcmc_sampler):
#        print(i)
        if (i+1)%500 == 0:
#            print('making sampler')
            temp_sampler = copy.copy(sampler)
            del temp_sampler.__dict__['lnprobfn']
            del temp_sampler.__dict__['pool']

            OutputData = {}
            OutputData['prefilename']=PreWalkingPickleFilename
            OutputData['sampler'] = temp_sampler
            pickle.dump(OutputData, open(FitOutput, 'wb'))

            del temp_sampler
        pass
#    for i, l_iterator_values in enumerate(mcmc_sampler):
#        print(i)
#        if (i+1)%500:
#            temp_sampler = sampler
#            del temp_sampler.__dict__['lnprobfn']
#            del temp_sampler.__dict__['pool']
#
#            OutputData = {}
#            OutputData['prefilename']=PreWalkingPickleFilename
#            OutputData['sampler'] = temp_sampler
#            pickle.dump(OutputData, open(FitOutput, 'wb'))

#sampler.run_mcmc(p0, niterations)



# print(sampler.chain[:,:,:]) 
# the sampler.chain is just a list of each (walker, iterator, dim) value


# delete un-pickleable parts of sampler: lnprobfn and pool
del sampler.__dict__['lnprobfn']
del sampler.__dict__['pool']

OutputData = {}
OutputData['prefilename']=PreWalkingPickleFilename
#OutputData['ndim'] = ndim
#OutputData['nwalkers'] = nwalkers
#OutputData['niterations'] = niterations
#OutputData['chain'] = sampler.chain
#OutputData['acceptance_fraction'] = sampler.acceptance_fraction
OutputData['sampler'] = sampler
pickle.dump(OutputData, open(FitOutput, 'wb'))

print(str((time.time() - StartingTimeFit)/3600.))
