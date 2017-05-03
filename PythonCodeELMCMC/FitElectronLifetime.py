################################
### Will use emcee for fitting
################################
import ElectronLifetimeTrend
from ElectronLifetimeTrend import *

import MyHistorianLib
from MyHistorianLib import GetUnixTimeFromTimeStamp

import FormPars
from FormPars import FormPars

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

# setting the parameters
MinUnixTime = GetUnixTimeFromTimeStamp(FormPars.MinTimeStamp)
MaxUnixTime = GetUnixTimeFromTimeStamp(FormPars.MaxTimeStamp)
default_pars = FormPars.GetDefaultPars()
#default_pars = [
#             3.7e-3, # attaching rate from literature
#             7.12997918e+03, # initial GXe concentration
#             5.09257628e+01, # initial LXe concentration
#             4.02212420e-01, # impurity attaching prob for vaporization
#             4.10390827e-01, # impurity attaching prob for condensation
#             1.90294571e+02, # GXe volume outgassing, in unit of kg/day
#             1.94252374e+02, # LXe volume outgassing, in unit of kg/day
#             [1465937520, 1468597800, 1479772379], # time for the impurity change, after correction
#             [0, 0, 0],
#             [1.00613537e-04, 3.36014333e-05, 1.0e-6],
#             [[1471880000, 1472800000, -100.]],
#             1000., # GXe outgassing linear decreasing constant, in days.
#             1000., # LXe outgassing linear decreasing constant, in days.
#             [[1487670000, 1000.], ], # additional LXe outgassing linear decreasing 
#             [1480317149, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#             [1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#             [1485351200, 1486272625, 0.9], # periods when getter is suspected to have lowered efficiency after earthquake
#            ]

# parameter selection
# select which parameters to fit
# 160609: 5 parameters
# 160809: 8 parameters
# 160812: 9 parameters
# 160902: 8 parameters
# 160927: 1 parameter
# 161003: 9 parameters
# 161021: 11 parameters
# 161212: 12 parameters with one more for the power glitch @ 11-21
# 161212: 13 parameters with an deficiency for getter after 11-28
# 170207: 14 parameters
# 170213: 15 parameters
# 170331: 14 parameters
# 170402: 16 parameters
#def FormPars(x):
#    if len(x)<16:
#        return default_pars
#    print("x=")
#    print(x)
#    IfOutOfBoundary = False
#    for i, y in enumerate(x):
#        if y<0 and (not i==9):
#            IfOutOfBoundary = True
#        if i==9 and y>0:
#            IfOutOfBoundary = True
#        if (i==2 or i==3 or i==13 or i==14 or i==15) and y>1:
#            IfOutOfBoundary = True
#    pars = default_pars
#    pars[1] = x[0] # initial GXe concentration
#    pars[2] = x[1] # initial LXe concentration
#    pars[3] = x[2] # vaporization attaching prob
#    pars[4] = x[3] # condensation attaching prob
#    pars[5] = x[4] # GXe outgassing
#    pars[6] = x[5] # LXe outgassing
#    pars[9][0] = x[6] # the amount impurity changed during power event
#    pars[9][1] = x[7] # the amount impurity changed during LN2 test @ July 15th
#    pars[9][2] = x[8] # the amount impurity changed during power glitch @ Nov. 21th
#    pars[10][0][2] = x[9] # the amount of outgassing in GXe changing due to gas-only flow
#    pars[11] = x[10] # GXe outgassing exponential decreasing constant, in days.
#    pars[12] = x[11] # LXe outgassing exponential decreasing constant, in days.
#    pars[13][0][1] = x[12] # LXe outgassing linear decreasing constant
#    pars[14][2] = x[13] # lowered efficiency
#    pars[15][2] = x[14] # lowered efficiency for Rn calibration during Christmas
#    pars[16][2] = x[15] # lowered efficiency for after earthquake in January
#    return (pars, IfOutOfBoundary)

# initial parameters

x0, x0_steps = FormPars.GetInitialParametersMCMC()
#x0 = np.array([5.82212635e+03, 
#                       6.55483847e+01,
#                       7.04240787e-01,
#                       2.20849318e-01,
#                       4.09610831e+02,
#                       6.34534001e+01,
#                       1.00e-5,
#                       2.77e-6,
#                       2.77e-7,
#                       -90,
#                       1000.,
#                       1000.,
#                       1000.,
#                       0.95,
#                       0.2,
#                       0.9,
#                      ])
#x0_steps = np.array([
#                                 5e2,
#                                 6.,
#                                 0.07,
#                                 0.02,
#                                 4.,
#                                 0.6,
#                                 1.e-6,
#                                 2.77e-7,
#                                 2.77e-8,
#                                 4.5,
#                                 1.,
#                                 1.,
#                                 1.,
#                                 0.05,
#                                 0.02,
#                                 0.05,
#                                ])

                                      

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
    contents = line2[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    #if was calculated using pax_v6.2.0 or younger
    if unixtime < 1478000000 or (unixtime > 1484900000 and unixtime < 1486100000):
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
    pars, IfOutOfBoundary = FormPars(x)
    if IfOutOfBoundary:
        print("Out of boundary!")
        return -np.inf
    global pElectronLifetimeTrend
    pElectronLifetimeTrend.SetParameters(pars)
    print("==== The parameters: =====")
    print(pElectronLifetimeTrend.GetParameters())
    LnL = LnLikeData()
    print(" LnL = "+str(LnL))
    print("=====================")
    print("Takes: "+str(time.time()-start_time) + " sec")
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
        p0.append([np.random.normal(a,b) for a, b in zip(x0, x0_steps)])

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, LnLike, threads=4)
sampler.run_mcmc(p0, niterations)

# print(sampler.chain[:,:,:]) 
# the sampler.chain is just a list of each (walker, iterator, dim) value

import pickle

OutputData = {}
OutputData['prefilename']=PreWalkingPickleFilename
OutputData['ndim'] = ndim
OutputData['nwalkers'] = nwalkers
OutputData['niterations'] = niterations
OutputData['chain'] = sampler.chain
OutputData['acceptance_fraction'] = sampler.acceptance_fraction
pickle.dump(OutputData, open(FitOutput, 'wb'))

print(str((time.time() - StartingTimeFit)/3600.))
