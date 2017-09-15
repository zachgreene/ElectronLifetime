import sys
import pickle
import corner
import numpy as np
from tqdm import tqdm
from MCMC_Tools import PrintParameterQuantilesWiki

import matplotlib.pyplot as plt

import FormPars

if len(sys.argv) < 2:
    print('========== Syntax ===========')
    print('python PlotCornerAndGelmanRubin.py')
    print('<mcmc pickle file>')
    print('<burn in cutoff (optional, default=200>')
    exit()

BurnInCutOff = 200

sPathToPickleMCMC = sys.argv[1]

if len(sys.argv) > 2:
    BurnInCutOff = int(sys.argv[2])



MCMC_Results = pickle.load(open(sPathToPickleMCMC, 'rb'))
sampler = MCMC_Results['sampler']
ndim = sampler.__dict__['dim']
niterations = sampler.__dict__['iterations']
nwalkers = sampler.__dict__['k']
chain = sampler.__dict__['_chain']

TotalEvents = chain.shape[1]

ParInfo = FormPars.GetParInfo()
ParNames = [ParInfo[i][0] for i in range(ndim)] # parameter names
ParDescs = [ParInfo[i][1] for i in range(ndim)] # parameter description
ParUnits = [ParInfo[i][2] for i in range(ndim)] # parameter units

#print(chain.shape)

assert TotalEvents == niterations
assert len(ParNames) == chain.shape[2]



#####################################
########  make corner plot  #########
#####################################

samples = chain[:, -BurnInCutOff:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=ParNames, quantiles=[0.16, 0.5, 0.84], show_titles=True)




#####################################
########  Gelman-Rubin Statistic  ###
#####################################

BatchSize = int(TotalEvents/40)
NumBatches = int(TotalEvents/BatchSize/2)
aForPlotting = [2*i*BatchSize for i in range(NumBatches)]

# define Gelman-Rubin for each variable
dGelmanRubin = {}
for Par in ParNames:
    dGelmanRubin[Par] = [0 for i in range(NumBatches)]


for i in tqdm(range(ndim)):
    Par = ParNames[i]
    for j in range(1, NumBatches + 1):
        NumEventsInBatch = j * BatchSize

        aMeans = np.mean(chain[:, j*BatchSize:2*j*BatchSize, i], axis=1)
        aVars = np.var(chain[:, j*BatchSize:2*j*BatchSize, i], axis=1, ddof=1)
    
        MeanOfMeans = np.mean(aMeans)
    
        B = NumEventsInBatch / (nwalkers - 1) * np.sum((aMeans - MeanOfMeans)**2) # inter-chain variance
#        W = 1. / nwalkers / (NumEventsInBatch - 1) * np.sum(aVars)          # intra-chain variance
        W = 1. / nwalkers * np.sum(aVars)          # intra-chain variance
    
        SigmaSquared = (NumEventsInBatch - 1) / (NumEventsInBatch) * W + B / NumEventsInBatch #weighted average variance
        V = SigmaSquared + B / (nwalkers * NumEventsInBatch)
    
        dGelmanRubin[Par][j-1] = (V / W)**0.5

PrintParameterQuantilesWiki(samples, [16, 50, 84], ParNames, ParDescs, ParUnits)

fig = plt.figure(figsize=(15, 20))
colors = plt.get_cmap('Paired')(np.linspace(0, 1., ndim))

for i, Par in enumerate(ParNames):
    plt.plot(aForPlotting, dGelmanRubin[Par], label=Par, color=colors[i], linewidth=3)

plt.hlines([1.1], 0, TotalEvents, color='k', linestyle='--')
plt.xlabel('Iterations')
plt.ylabel('Gelman-Rubin Statistic')
plt.gca().set_yscale('log')
plt.legend()
plt.show()
