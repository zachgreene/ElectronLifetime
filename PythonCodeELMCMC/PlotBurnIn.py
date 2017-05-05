import pickle
import matplotlib.pyplot as plt
import sys, os
import numpy as np

import FormPars

if len(sys.argv) < 2:
    print('========== Syntax ===========')
    print('python PlotBurnIn.py')
    print('<mcmc pickle file>')
    print('<path to figure directory (optional)>')
    exit()


sPathToPickleMCMC = os.path.abspath(sys.argv[1])

if len(sys.argv) > 2:
    sPathToSaveDir = sys.argv[2]
else:
    sPathToSaveDir = '/'.join(sPathToPickleMCMC.split('/')[:-2]) + '/Figures/'

MCMC_Name = sPathToPickleMCMC.split('/')[-1].split('.')[0]

############# load MCMC sampler #############
MCMC_Results = pickle.load(open(sPathToPickleMCMC, 'rb'))
sampler = MCMC_Results['sampler']
ndim = sampler.__dict__['dim']
niterations = sampler.__dict__['iterations']
nwalkers = sampler.__dict__['k']
chain = sampler.__dict__['_chain']

TotalEvents = chain.shape[1]


############# get parameter information #############
ParInfo = FormPars.GetParInfo()
ParNames = [ParInfo[i][0] for i in range(ndim)] # parameter names
ParDescs = [ParInfo[i][1] for i in range(ndim)] # parameter description
ParUnits = [ParInfo[i][2] for i in range(ndim)] # parameter units



############# set plot information #############
colors = ['k', 'b','r','g','y','darkblue', 'darkgreen', 'chocolate', 'm', 'gold', 'purple', 'c', 'violet', 'r', 'k', 'darkorange', 'forestgreen', 'lightseagreen', 'rebeccapurple']

Iterations = np.linspace(1, niterations, niterations)

fig, ax = plt.subplots(int(ndim/2. + 0.5), 2, sharex=True, figsize=(20,20))
ax = ax.T.reshape((-1))


############# plot walkers vs. iteration for each parameter #############
for i in range(ndim):
    ax[i].plot(Iterations, chain[:, :, i].T, color=colors[i], linewidth=1, alpha=0.1)
    # set units for ylabel if exist
    if ParUnits[i] != '':
        yLabel = ParNames[i] + ' [' + ParUnits[i] + ']'
    else:
        yLabel = ParNames[i]
    ax[i].set_ylabel(yLabel)
    # label x-axis for bottom plots
    if (i + 1) == len(ax)/2:
        ax[i].set_xlabel('Iterations')
    ax[i].yaxis.get_major_formatter().set_powerlimits((-3,3))

ax[-1].set_xlabel('Iterations')

plt.savefig(sPathToSaveDir + MCMC_Name + '_BurnIn.png', format='png')
plt.show()
