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

print('\nFitting Electron Lifetime between ' + FormPars.GetMinTimeStamp() + ' and ' + FormPars.GetMaxTimeStamp() + '\n')

# setting the parameters
MinUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMinTimeStamp())
MaxUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMaxTimeStamp())
default_pars = FormPars.GetDefaultPars()

# initial parameters
x0, x0_steps = FormPars.GetInitialParametersMCMC()

# The main Light yield Trend
pElectronLifetimeTrend = MyElectronLifetimeTrend(HistorianFile, MinUnixTime, MaxUnixTime, default_pars)

# pre-walking
nwalkers = int(sys.argv[5])
niterations = int(sys.argv[6])
PreWalkingPickleFilename = "NoneExist"
if len(sys.argv)>7:
    PreWalkingPickleFilename = sys.argv[7]

ElectronLifetimeData = {}
#############################
# fill in the data
#############################

# load eelectron lifetime data

UnixTimes, UnixTimeErrors, Values, ValueErrors = LoadFitData('SingleScatter', PathToFile=ElectronLifetimeDataFile)
RnUnixtimes, RnUnixtimeErrors, RnELifeValues, RnELifeValueErrors = LoadFitData('Rn', PathToFile=RnElectronLifetimeDataFile)

## Load the Rn data

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
