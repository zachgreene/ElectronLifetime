import root_pandas
import hax
import pandas as pd
import os, glob
import sys
import datetime, time
from pymongo import MongoClient
import socket
from Cuts import *
from Tools import *
import BasicPlots
import MongoDB_Query
import hax

import matplotlib.pyplot as plt
from matplotlib import colors
# from matplotlib.colors import LogNorm

import multihist
import numpy as np
from scipy.optimize import curve_fit

sPathToFigureSaveDirectory = ''
bShowBasicFigs = 0
pax_version = '6.4.2'
#sPathToFilelist = '/home/zgreene/xenon1t/ElectronLifetime/Filelists/170109_0732_24hrs.txt'

if len(sys.argv) < 2:
    print('====== Usage ======')
    print('python ExtractElectronLifetime.py')
    print('<list>')
    print('<lifetime output txt file>')
    print('<figure directory>')
    print('<show basic figs (0 for no, 1 for yes)>')
    print('<pax_version>')
    exit()


sPathToFilelist = sys.argv[1]
FilelistName = sPathToFilelist.split('/')[-1].split('.')[0]
sPathToLifetimesFile = sys.argv[2]
if len(sys.argv) > 3:
    sPathToFigureSaveDirectory = sys.argv[3]
    sPathToFigureSaveDirectory += FilelistName
    if not os.path.isdir(sPathToFigureSaveDirectory):
        os.mkdir(sPathToFigureSaveDirectory)
    sPathToFigureSaveDirectory += '/'
if len(sys.argv) > 4:
    bShowBasicFigs = int(sys.argv[4])
if len(sys.argv) > 5:
    pax_version = sys.argv[5]



#sPathToDatasetList = '/home/zgreene/xenon1t/ElectronLifetime/ELifeTXTs/ElectronLifetimeDatasets_pax_v'+pax_version+'.txt'
sPathToDatasetList = '/home/zgreene/xenon1t/ElectronLifetime/ELifeTXTs/ElectronLifetimeDatasets.txt'

aPathToMinitrees = ['/project/lgrandi/xenon1t/minitrees/pax_v' + pax_version + '/']
aPathToMinitrees.append('/project2/lgrandi/xenon1t/minitrees/pax_v' + pax_version + '/')
aPathToMinitrees.append('/project/lgrandi/zgreene/xenon1t/minitrees/pax_v' + pax_version + '/')

# percentiles to find lifetime
qs = [16,50,84]



############################################
###### read in filelist and load data ######
############################################

print('Loading files from ' + FilelistName)
RunInfo = ReadFilelist(sPathToFilelist)
data = LoadData(RunInfo.name, aPathToMinitrees)

data = data.merge(RunInfo, left_on='run_number', right_on='number', how='left')

UniqueSources = data.source.unique()
print('Sources present during runs: ' + ', '.join(UniqueSources))



############################################
################# get cuts #################
############################################

CutLimits = GetInitialCuts(pax_version, RunInfo.source)
#Xs1 = EvalCut(data, CutLimits, 'Xs1')
#Xs2 = EvalCut(data, CutLimits, 'Xs2')
#Xs1asym = EvalCut(data, CutLimits, 'Xs1asym')
#Xs2asym = EvalCut(data, CutLimits, 'Xs2asym')
#Xr = EvalCut(data, CutLimits, 'Xr')
#Xz = EvalCut(data, CutLimits, 'Xz')
CutLimits['MaxValueDt'] = MaxValueDt = GetMaxDriftTime(data[np.isfinite(data.dt)].dt)
#CutLimits['MaxValueDt'] = MaxValueDt = 727.5
print(CutLimits['MaxValueDt'])
dataOrig = data
data = hax.cuts.range_selection(data, 's1', (CutLimits['LowLimitS1'], CutLimits['UpLimitS1']))
data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))
data = hax.cuts.range_selection(data, 's1_aft', (CutLimits['LowLimitAsymS1'], CutLimits['UpLimitAsymS1']))
data = hax.cuts.range_selection(data, 's2_aft', (CutLimits['LowLimitAsymS2'], CutLimits['UpLimitAsymS2']))
data = hax.cuts.selection(data, data.dt < MaxValueDt, 'dt < %i' %MaxValueDt)
data = hax.cuts.selection(data, data.r < CutLimits['RadialLimit'], 'r2 < %i' %(CutLimits['RadialLimit'])**2)
#print(CutLimits['MaxValueDt'])

#Xcuts = Xs1 & Xs2 & Xs1asym & Xs2asym & Xr & Xz




############################################
###### generate basic plots with cuts ######
############################################
figs = []

# make plots if either will be saved or shown, otherwise not
if sPathToFigureSaveDirectory != '' or bShowBasicFigs:
    figs.append(BasicPlots.PlotS1S1Asym(dataOrig, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotS1S2(dataOrig, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotS1Asym(dataOrig, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotS2Asym(dataOrig, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotRadius(dataOrig, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotRadiusFlat(data, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotDepth(dataOrig, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotDepthFlat(data, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    figs.append(BasicPlots.PlotS1S1AsymWithCuts(data, CutLimits, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory))
    
#    if not bShowBasicFigs:
#        for fig in figs:
#            plt.close(fig)
   
#    plt.show()


############################################
########## fit electon lifetime  ###########
############################################
MaxValueDt = int(CutLimits['MaxValueDt'])
BinsOmitLow = CutLimits['BinsOmitLow']
BinsOmitUp = CutLimits['BinsOmitUp']

# exclude points on either side of this factor multiplied by prelim fit
FactorBelow = CutLimits['FactorBelow']
FactorAbove = CutLimits['FactorAbove']

TimePoints = np.linspace(0, MaxValueDt, MaxValueDt + 1)
x = data.dt
y = data.s2
xBins = np.linspace(0, MaxValueDt, 100)
yBins = np.linspace(0, 150000, 100)

# get preliminary fit
try:
    PrelimOpts = GetPreliminaryFit(x, y, xBins, qs=qs, FitType='Exp', BinsOmitLow=BinsOmitLow, BinsOmitUp=BinsOmitUp) 
except:
    print('error with prelim 100 bins, trying with 80')
    PrelimOpts = GetPreliminaryFit(x, y, np.linspace(0, MaxValueDt, 80), qs=qs, FitType='Exp', BinsOmitLow=BinsOmitLow, BinsOmitUp=BinsOmitUp) 
    print('prelim fit worked with 80 bins')

#Xs2dtLow = SetS2DtLow(data, PrelimOpts, FactorBelow)
#Xs2dtUp = SetS2DtUp(data, PrelimOpts, FactorAbove)
#Xs2dt = Xs2dtLow & Xs2dtUp

p0, p1 = PrelimOpts
dataOmit = hax.cuts.selection(data, (data.s2 < FactorBelow * p0 * np.exp(-data.dt / p1)) |
                                    (data.s2 > FactorAbove * p0 * np.exp(-data.dt / p1)),
                                    's2-dt excluded')
data = hax.cuts.selection(data, data.s2 > FactorBelow * p0 * np.exp(-data.dt / p1), 's2-dt lower bound')
data = hax.cuts.selection(data, data.s2 < FactorAbove * p0 * np.exp(-data.dt / p1), 's2-dt upper bound')

# draw s1 s2 plot
figS1S2 = BasicPlots.PlotS1S2WithCuts(data, sPathToFigureSaveDirectory=sPathToFigureSaveDirectory)
#plt.close(figS1S2)

# to draw as s2-dt boundaries on plot
PrelimFitLine = FitExponential(TimePoints, *PrelimOpts)

### perform final fit ###
#try:
##    BinCenters, S2Percentiles = GetBinPercentiles(data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['dt'],
##                                             data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['s2'],
##                                             xBins=np.linspace(0,MaxValueDt,100), qs=qs)
#    BinCenters, Medians, MADs = GetBinMediansAndMADs(data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['dt'],
#                                             data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['s2'],
#                                             xBins=np.linspace(0,MaxValueDt,100))
#except:
#    print('error with 100 bins, trying with 80')
##    BinCenters, S2Percentiles = GetBinPercentiles(data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['dt'],
##                                             data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['s2'],
##                                             xBins=np.linspace(0,MaxValueDt,80), qs=qs)
#    BinCenters, Medians, MADs = GetBinMediansAndMADs(data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['dt'],
#                                             data[Xs1 & Xs1asym & Xs2asym & Xr & Xz & Xs2dt]['s2'],
#                                             xBins=np.linspace(0,MaxValueDt,100))
#    print('worked with 80 bins')


sPathToFileSave = sPathToFigureSaveDirectory + 'GaussFits_' + FilelistName + '.png'

BinCenters, BinErrs, OptParams, OptCovs = GetBinGaussianParameters(data.dt,
                                             data.s2,
                                             xBins=np.linspace(0,MaxValueDt,20),
                                             sPathToFileSave=sPathToFileSave,
                                             bUseROOT=True)


#LowErrs, Medians, UpErrs = S2Percentiles.T
#AvgErrs = (UpErrs - LowErrs)/2.

#LowErrs = Medians - MADs
#UpErrs = Medians + MADs
#AvgErrs = MADs

Means = OptParams.T[1]
AvgErrs = np.sqrt(OptCovs.T[1][1])
UpErrs = Means + AvgErrs
LowErrs = Means - AvgErrs
#popt, pcov = curve_fit(FitExponential, BinCenters[BinsOmitUp:-BinsOmitLow],
#                          Medians[BinsOmitUp:-BinsOmitLow],
#                          [1e5, 400], sigma=AvgErrs[BinsOmitUp:-BinsOmitLow])

#popt, pcov = curve_fit(FitExponential, BinCenters[BinsOmitUp:-BinsOmitLow],
#                        Means[BinsOmitUp:-BinsOmitLow],
#                        [1e5, 400], sigma=AvgErrs[BinsOmitUp:-BinsOmitLow])


bStatusFail, popt, pcov = FitExponentialROOT(BinCenters[BinsOmitUp:-BinsOmitLow],
                        Means[BinsOmitUp:-BinsOmitLow],
                        BinErrs[BinsOmitUp:-BinsOmitLow],
                        AvgErrs[BinsOmitUp:-BinsOmitLow],
                        [1e5, 400])

Lifetime = popt[1]
LifetimeErr = (pcov[1][1])**0.5

# get time and error
StartTimeUnix, EndTimeUnix = MongoDB_Query.GetStartStopTimes(RunInfo.name)
MeanTime = 0.5*(StartTimeUnix + EndTimeUnix)
ErrTime = 0.5*(EndTimeUnix - StartTimeUnix)



############################################
########## draw electon lifetime  ##########
############################################
#figLifetime = plt.figure(figsize=(15,8))

#if not LifetimeLogPlot:
#    # draw points excluded by s2-dt boundaries
#    plt.hist2d(data[Xcuts & (Xs2dt==0)]['dt'],
#                   data[Xcuts & (Xs2dt==0)]['s2'],
#                   bins=[np.linspace(0,MaxValueDt,100), np.linspace(0,150000,100)],
#                   norm=colors.LogNorm(), cmap='viridis', alpha=0.5)
#    
#    # draw points that will be used in final fit
#    plt.hist2d(data[Xcuts & Xs2dt]['dt'],
#                   data[Xcuts & Xs2dt]['s2'],
#                   bins=[np.linspace(0,MaxValueDt,100), np.linspace(0,150000,100)],
#                   norm=colors.LogNorm(), cmap='viridis')

figLifetime = GetLifetimeLogPlot(data, dataOmit, MaxValueDt, S2Min=5000, S2Max=200000)
# draw s2-dt cut lines
plt.plot(TimePoints, FactorAbove*PrelimFitLine,
            linestyle='dashed', color='purple', label='Remove events above distribution')
plt.plot(TimePoints, FactorBelow*PrelimFitLine,
            linestyle='dashed', color='blue', label='Remove events below distribution')

# draw dt cuts from bins
plt.vlines([np.linspace(0,MaxValueDt,100)[BinsOmitUp], np.linspace(0,MaxValueDt,100)[-BinsOmitLow]],
           plt.ylim()[0], plt.ylim()[1], color='turquoise', linestyle='dashed', label='Lifetime Fit Cuts')


# draw percentiles
plt.errorbar(BinCenters[BinsOmitUp:-BinsOmitLow],
                        Means[BinsOmitUp:-BinsOmitLow],
                        xerr=BinErrs[BinsOmitUp:-BinsOmitLow],
                        yerr=AvgErrs[BinsOmitUp:-BinsOmitLow],
                        fmt='o', color='orangered', capsize=0,
                        markeredgewidth=0, linewidth=2)
#plt.step(BinCenters, Means, color='lime', linestyle='-')
#plt.step(BinCenters, LowErrs, color='darkturquoise', linestyle='-')
#plt.step(BinCenters, UpErrs, color='darkturquoise', linestyle='-')
# draw best fit lifetime
plt.step(TimePoints, FitExponential(TimePoints, *popt), linewidth=2, color='r', linestyle='-')

plt.xlabel('dt $[\mu s]$')
plt.ylabel('S2 [PE]')
#plt.title(FilelistName)
plt.legend()
plt.colorbar(label='Counts')
plt.tight_layout()

plt.gca().set_yscale('log')
plt.xlim(0, MaxValueDt)
plt.ylim(15000, 200000)
LifetimeText = 'Electron Lifetime = %.2f $\pm$ %.2f $\mu s$' % (Lifetime, LifetimeErr)
#plt.text(40, 0.9*plt.ylim()[1], LifetimeText, color='r', fontsize=18, weight='bold', alpha=1)
#plt.text(40, 0.84*plt.ylim()[1], '@ %.3f $\pm$ %.3f' % (MeanTime, ErrTime), color='r', fontsize=18, weight='bold', alpha=1)
plt.text(40, 1.4*plt.ylim()[0], LifetimeText, color='r', fontsize=18, weight='bold', alpha=1)
plt.text(40, 1.25*plt.ylim()[0], '@ %.3f $\pm$ %.3f' % (MeanTime, ErrTime), color='r', fontsize=18, weight='bold', alpha=1)

# if electron lifeitme fit did not fail
if not bStatusFail:
    # keep list of all runs used to compute electron lifetime
    InsertRunsToDatasetList(RunInfo.name, sPathToDatasetList)
    ## insert electron lifetime to file
    ValuesToInsert = (MeanTime, ErrTime, Lifetime, LifetimeErr)
    InsertLifetimeToFile(ValuesToInsert, sPathToLifetimesFile, bOverwrite=False)

print('Electron Lifetime = %.2f +/- %.2f' % (Lifetime, LifetimeErr))

#print(pax_version, LifetimeText)
#if bSaveFigs:
plt.savefig(sPathToFigureSaveDirectory + 'LifetimeFit_' + FilelistName + '_Gauss.png')
#plt.close(figLifetime)
plt.show()
