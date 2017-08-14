import ElectronLifetimeTrend
from ElectronLifetimeTrend import *

#import MyHistorianLib
#from MyHistorianLib import GetUnixTimeFromTimeStamp

import Tools
from Tools import *

import FormPars

import numpy as np
from scipy.optimize import minimize
from numpy.linalg import inv

if len(sys.argv)<2:
    print("======= Syntax ========")
    print("python PredictElectronLifetime.py")
    print("<Historian file>")
    print("<Fit result pickle>")
    print("<Prediction txt output>")
    print("<burn-in iteraction cut>")
    print("NOTICE: do remember to add \"-r\" at the end")
    exit()

HistorianFile = sys.argv[1]
FitResultInput = sys.argv[2]
PredictionOutputFile = sys.argv[3]
BurnInCutOff = int(sys.argv[4])


NumOfInterpolation = 5000
NumOfTrials = 1000


# setting the parameters
MinUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMinTimeStamp())
MaxUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMaxTimeStamp()) + 50*24*3600
default_pars = FormPars.GetDefaultPars()


#############################
## Get the MCMC result
#############################
import pickle
MCMCResults = pickle.load(open(FitResultInput, 'rb'))
sampler = MCMCResults['sampler']
ndim = sampler.__dict__['dim']
niterations = sampler.__dict__['iterations']
nwalkers = sampler.__dict__['k']
print(ndim, nwalkers, niterations)

# Cutting Burn-In and unreasonable region
# the value shall already be from "PlotElectronLifetime.py"
#samples_cut = GetBurnInCutoffSamples(sampler.__dict__['_chain'], int(BurnInCutOff))
samples_cut = GetBurnInCutoffSamplesV2(sampler.__dict__['_chain'], int(BurnInCutOff), int(BurnInCutOff)+199)

#####################################
## Calculate the best fit from MCMC results
#####################################
mean = np.average(samples_cut, axis=0)
print(mean)


Pars, IfSth=FormPars.FormPars(mean)
#MaxUnixTime = LastPointUnixTime + 60.*3600.*24. # 2 month after the last data point
# The main Light yield Trend
pElectronLifetimeTrend = MyElectronLifetimeTrend(HistorianFile, MinUnixTime, MaxUnixTime, Pars)

# get the graphs 
UnixTimes = np.linspace(MinUnixTime, MaxUnixTime, NumOfInterpolation)
# include times of discontinuity
ImpactfulUnixtimes = FormPars.GetImpactfulUnixtimes()
for ImpactfulUnixtime in ImpactfulUnixtimes:
	if ImpactfulUnixtime-1 not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime-1]), ImpactfulUnixtime-1)
	if ImpactfulUnixtime not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime]), ImpactfulUnixtime)
	if ImpactfulUnixtime+1 not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime+1]), ImpactfulUnixtime+1)

#############################
# calculate:
# 1) Maximum
# 2) Reaching date for 90% capacity
# 3) Reaching date for 500us
#############################
# The electron lifetime Trend for varied parameters

# function for checking the list percentile with np.inf in it
def PercentileWithInf(Values, deviation):
    TotalCounter = len(Values)
    InfCounter = 0
    NewValues = []
    # first get the new list and counter
    for value in Values:
        if value==np.inf:
            InfCounter += 1
        else:
            NewValues.append(value)
    Fraction = (1. - float(InfCounter) / float(TotalCounter)) * 100.
    if deviation>Fraction:
        return np.inf
    deviation_new = deviation / Fraction *100.
    return np.percentile(NewValues, deviation_new, axis=0)

# Initial the 1 sigma lower/upper of the trend
Trends = []
for i in range(len(UnixTimes)):
    Trends.append([])


# pickup randomly the pars in the MCMC walker
Pars_Trials = PickupMCMCPosteriors(samples_cut, NumOfTrials)
for i, pars_random in enumerate(Pars_Trials):
    print("i = "+str(i))
    pars, IfSth = FormPars.FormPars(pars_random)
    pElectronLifetimeTrend.SetParameters(pars)
    for j, unixtime in enumerate(UnixTimes):
        Trends[j].append(pElectronLifetimeTrend.GetElectronLifetime(unixtime))
StandardDeviation1Sigma = [15.4, 50., 84.6]
Taus = []
LowerBoundaries = []
UpperBoundaries = []
for i, deviation in enumerate(StandardDeviation1Sigma):
    for Values in Trends:
        Boundary = PercentileWithInf(Values, deviation)
        if i==0:
            LowerBoundaries.append(Boundary)
        elif i==1:
            Taus.append(Boundary)
        else:
            UpperBoundaries.append(Boundary)
#print(MaximumELifes)




######### moved from bottom because terminal cannot show root plots ##########
fout = open(PredictionOutputFile+".txt", 'w')
for unixtime, tau, lower, upper in zip(UnixTimes, Taus, LowerBoundaries, UpperBoundaries):
    fout.write(str(unixtime)+"\t\t")
    fout.write(str(tau)+"\t\t")
    fout.write(str(lower)+"\t\t")
    fout.write(str(upper)+"\t\t")
    fout.write("\n")
fout.close()
#######################################################################
