import datetime
import time
import numpy as np
import pandas as pd
import os
import root_pandas
from pymongo import MongoClient
import socket
from scipy.optimize import curve_fit
import ROOT
import matplotlib.pyplot as plt
from matplotlib import colors

def GetUnixtime(Filename):
    contents = Filename.split("/")
    Filename = contents[-1]
    Year = "20"+Filename[0:2]
    Month = Filename[2:4]
    Day = Filename[4:6]
    Hour = Filename[7:9]
    Minute = Filename[9:11]
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), 0)
    return time.mktime(dt.timetuple())


def GetRunNamesFromMinitreeRunList(MinitreeRunList):
    RunNames = []
    for Minitree in MinitreeRunList:
        contents = Minitree.split('/')[-1].split('_')
        RunName = contents[0] + '_' + contents[1]
        if len(RunName) < 2:
            continue
        RunNames.append(RunName)

    return RunNames



def ReadTotalRunsList(sPathToFilelist):
    fin = open(sPathToFilelist, 'r')
    lines = fin.readlines()
    fin.close()

    RunNames = []
    for line in lines:
        if len(line) < 2:
            continue
        line = line[:-1]
        contents = line.split('\t\t')
        RunNames.append(contents[0])

    return RunNames



def ReadFilelist(sPathToFilelist):
    fin = open(sPathToFilelist, 'r')
    lines = fin.readlines()
    fin.close()

    RunNames, RunNumbers, RunSources, RunStartTimes, RunEndTimes = [], [], [], [], []
    for line in lines:
        if len(line) < 2:
            continue
        line = line[:-1]
        contents = line.split('\t\t')
        RunNames.append(contents[0])
        RunNumbers.append(int(contents[1]))
        RunSources.append(contents[2])
        RunStartTimes.append(float(contents[3]))
        RunEndTimes.append(float(contents[4]))

    dRunInfo = dict(
                    name = RunNames,
                    number = RunNumbers,
                    source = RunSources,
                    start = RunStartTimes,
                    end = RunEndTimes
                    )

    RunInfo = pd.DataFrame.from_dict(dRunInfo)
    RunInfo = RunInfo[['name', 'number', 'source', 'start', 'end']]

#    return (RunNames, RunNumbers, RunSources)
    return RunInfo


def OutputFilelist(FilesInList, sPathToDatasetsList):
    fout = open(sPathToDatasetsList, 'w')

    for File in FilesInList:
        fout.write(File + '\n')

    fout.close()
    return


def LoadData(RunNamesToLoad, aPathToMinitrees, MinitreeType='Basics', ColumnsToLoad='*', Preselection='event_number != -1', verbose=False):

    FilesNotFound = []
    FileDfs = []       # each root file gets own df, then concat

    TotalEvents = 0

    for RunName in RunNamesToLoad:
        bFoundFile = False
        for sPathToMinitrees in aPathToMinitrees:
            file = sPathToMinitrees + RunName + '_' + MinitreeType + '.root'

        #   check if file exists
            if os.path.isfile(file):
                if verbose:
                    print('Loading ' + file)
                FileDf = root_pandas.read_root(file, columns=ColumnsToLoad)
                TotalEvents += len(FileDf)
                try:
                    FileDf = FileDf[FileDf.eval(Preselection)]
                    FileDfs.append(FileDf)
                    bFoundFile = True
                except:
                    pass
                break

#         if file does not exist
        if not bFoundFile:
            FilesNotFound.append(RunName)
            print('!!!! Unable to load ' + RunName + ' !!!!')

        if len(FileDfs) > 0:
            data = pd.concat(FileDfs)

    print('====> %i/%i root files found and loaded with %i events' % (len(FileDfs), len(RunNamesToLoad), len(data)))
    if Preselection != 'event_number != -1':
        PercentPassed = 100.*len(data)/TotalEvents
        print('====> %i/%i (%.2f%%) passed preselection cut \'%s\'' %(len(data), TotalEvents, PercentPassed, Preselection))

        
    if MinitreeType=='Basics':
        data['dt'] = data['drift_time']/1000
        data['r'] = np.power(np.power(data['x'], 2) + np.power(data['y'], 2), 0.5)
        data['s2_bottom'] = data['s2']*(1 - data['s2_area_fraction_top'])
        data['cs2_bottom'] = data['cs2']*(1 - data['s2_area_fraction_top'])
        data['s2aft_diff_up'] = data.s2_area_fraction_top-(0.68 + 0.02 / 650 * data.dt) 
        data['s2aft_diff_low'] = data.s2_area_fraction_top-((0.35 + 0.18* (1.0/(1+np.exp(-(data.dt-5)/20))))*(1+0.06/600*data.dt)) 
        data['s2s_diff'] = data.largest_other_s2 - (0.01*data.s2+40)
        data = data.rename(columns={'s1_area_fraction_top': 's1_aft', 's2_area_fraction_top': 's2_aft'})

    return data



def GetBinCenters(Bins):
    return 0.5*(Bins[1:] + Bins[:-1])


def GetMaxDriftTime(DriftTimes, bins=np.linspace(300,1000,701)):
    Counts, DriftTimeEdges = np.histogram(DriftTimes, bins=bins)
    DriftTimeCenters = GetBinCenters(DriftTimeEdges)
    # largest drift time should be max (outside of very beginning, so do not bin first 100 us)
    MaxDriftTime = DriftTimeCenters[Counts == np.max(Counts)][-1]
    return MaxDriftTime


def FitExponential(x, p0, p1):
    return p0*np.exp(-x/p1)


def FitExponential(x, p0, p1, p2=0, x0=0, bIncludeConst=False):
    if bIncludeConst:
        return p0*np.exp(-(x-x0)/p1) + p2
    else:
        return p0*np.exp(-x/p1) 




def FitExponentialROOT(xCenters, yValues, xValErrs, yValErrs, Pars, x0=0,  bIncludeConst=False):
    bStatusFail = False

    xMin = xCenters[0] - xValErrs[0]
    xMax = xCenters[-1] + xValErrs[-1]
    g = ROOT.TGraphErrors()

    for i in range(len(yValues)):
        g.SetPoint(i, xCenters[i], yValues[i])
        g.SetPointError(i, xValErrs[i], yValErrs[i])

    # choose which fit to perform
    if bIncludeConst:
        f = ROOT.TF1('f', '[0]*exp(-(x - %.f)/[1]) + [2]' %x0, xMin, xMax)
    else:
        f = ROOT.TF1('f', '[0]*exp(-x/[1])', xMin, xMax)

    f.SetParameters(*Pars)

    if Pars[1] > 0:
        f.SetParLimits(0, Pars[0] / 10., Pars[0] * 10.)
        f.SetParLimits(1, Pars[1] / 10., Pars[1] * 10.)
    elif Pars[1] < 0:
        f.SetParLimits(0, Pars[0] * 10., Pars[0] / 10.)
        f.SetParLimits(1, Pars[1] * 10., Pars[1] / 10.)

    ptr = g.Fit('f', 'RMELSN')
    if ptr.Status() % 10 != 0 or ptr.Status() % 100 != 0 or ptr.Status() % 1000 != 0:
        bStatusFail = True 

    print(ptr.Status(), type(ptr.Status()))

    # popt and pcov are different depending on if constant was included
    if bIncludeConst:
        popt = [f.GetParameter(0), f.GetParameter(1), f.GetParameter(2)]
        pcov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(9):
            pcov[int(i/3)][i%3] = ptr.GetCovarianceMatrix().GetMatrixArray()[i]

    else:
        popt = [f.GetParameter(0), f.GetParameter(1)]
        pcov = [[0, 0], [0, 0]]
        for i in range(4):
            pcov[int(i/2)][i%2] = ptr.GetCovarianceMatrix().GetMatrixArray()[i]

    print(popt, pcov)
    print(pcov[1][1]**0.5, f.GetParError(1))

    popt = np.asarray(popt)
    pcov = np.asarray(pcov)

    return (bStatusFail, popt, pcov)

def FitGaussian(x, p0, p1, p2):
    return p0*np.exp(-0.5*(x-p1)**2/p2**2)


def FitGaussianROOT(xEdges, yValues, Pars):
#    g = ROOT.TGraphErrors()
    xMin = min(xEdges)
    xMax = max(xEdges)

    if len(yValues) == 0:
        return (True, [0], [0])

    bStatusFail = False # whether or not fit failed
#    print(xMin, xMax)
    h = ROOT.TH1F('h', 'h', len(xEdges)-1, xMin, xMax)
#    xErr = (xCenters[1:] - xCenters[:-1])/2.
    for yValue in yValues:
#        print(yValue)
        h.Fill(yValue)
#        g.SetPoint(i, yCenters[i], Hist[i])
#        g.SetPointError(i, xErr[i], Hist[i]**0.5)

#    xMin = min(xCenters) - xErr
#    xMax = max(xCenters) + xErr

    # do preliminary fit first
    fPrelim = ROOT.TF1('fPrelim', 'gaus', xMin, xMax)
    fPrelim.SetParameters(*Pars)
    ptrPrelim = h.Fit('fPrelim', 'RMELSN')

    # follow with final fit
    xMinFinal = fPrelim.GetParameter(1) - 2*fPrelim.GetParameter(2)
    xMaxFinal = fPrelim.GetParameter(1) + 2*fPrelim.GetParameter(2)
    f = ROOT.TF1('f', 'gaus', xMinFinal, xMaxFinal)
    f.SetParameters(*Pars)

    ptr = h.Fit('f', 'RMELSN')

    print(ptr.Status(), type(ptr.Status()))
    if ptr.Status() % 10 != 0 or ptr.Status() % 100 != 0 or ptr.Status() % 1000 != 0:
        bStatusFail = True 

    popt = [f.GetParameter(0), f.GetParameter(1), f.GetParameter(2)]
    pcov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(9):
        pcov[int(i/3)][i%3] = ptr.GetCovarianceMatrix().GetMatrixArray()[i]

#    print(popt, pcov)
#    print(pcov[1][1]**0.5, f.GetParError(1))

    return (bStatusFail, popt, pcov)

    



def GetBinPercentiles(x, y, xBins, qs):
#   make sure is array, if not convert
    if isinstance(qs,int) or isinstance(qs,float):
        qs = [qs]

    xCenters = GetBinCenters(xBins)

    data = np.array([x, y])
    sortedData = [data[1][(data[0] > xBins[i]) & (data[0] <= xBins[i+1])] for i in range(len(xBins)-1)]
    Percentiles = np.asarray([np.percentile(sortedData[i], qs, axis=0) for i in range(len(sortedData))])

    return (xCenters, Percentiles)


def GetBinMediansAndMADs(x, y, xBins):
    xCenters = GetBinCenters(xBins)

    data = np.array([x, y])
    sortedData = [data[1][(data[0] > xBins[i]) & (data[0] <= xBins[i+1])] for i in range(len(xBins)-1)]
    Medians = np.asarray([np.median(sortedData[i], axis=0) for i in range(len(sortedData))])
    MADs = np.asarray([np.median(np.abs(sortedData[i] - Medians[i])) for i in range(len(sortedData))])

    return (xCenters, Medians, MADs)

def GetBinGaussianParameters(x, y, xBins, sPathToFileSave, bUseROOT=True, bWriteParsToFile=True):
    xCenters = GetBinCenters(xBins)
    xErr = (xCenters[1] - xCenters[0]) / 2.
    BinCenters = [] # will be equal to xCenters if bStatusFalse=False
    BinErrs = [] # will be equal to xErr if bStatusFalse=False

    data = np.array([x, y])
    sortedData = [data[1][(data[0] > xBins[i]) & (data[0] <= xBins[i+1])] for i in range(len(xBins)-1)]

    popts, pcovs = [], []
    plt.figure(figsize=(50,34))
    for iBin in range(len(sortedData)):
        yMin = np.min(sortedData[iBin])
        yMax = np.max(sortedData[iBin])
        Hist, yEdges = np.histogram(sortedData[iBin], np.linspace(yMin, yMax, 50))
#        print(yMin, yMax)
        print(Hist)
        yCenters = GetBinCenters(yEdges)
        if bUseROOT:
            bStatusFail, popt, pcov = FitGaussianROOT(yEdges, sortedData[iBin], Pars=[10,(yMax+yMin)/2,(yMax-yMin)/2])
        else:
            popt, pcov = curve_fit(FitGaussian, yCenters, Hist, [10,(yMax+yMin)/2,(yMax-yMin)/2])
            bStatusFail = False

        if not bStatusFail:
            BinCenters.append(xCenters[iBin])
            BinErrs.append(xErr)
            popts.append(popt)
            pcovs.append(pcov)

        # save figures of gauss
        plt.subplot(4,5,iBin+1)
        plt.errorbar(yCenters, Hist, yerr=Hist**0.5, fmt='o')
        aGaussS2s = np.linspace(yMin, yMax, 1000)
        plt.plot(aGaussS2s, FitGaussian(aGaussS2s, *popt))
        plt.text(1.2*plt.xlim()[0], 0.8*plt.ylim()[1], '$\mu_{S2}$ = %.1f $\pm$ %.1f pe' %(popt[1], pcov[1][1]**0.5),
                fontsize=30, weight='bold', color='r')

    plt.savefig(sPathToFileSave)
#    plt.close(plt.gcf())

    popts = np.asarray(popts)
    pcovs = np.asarray(pcovs)
    BinCenters = np.asarray(BinCenters)
    BinErrs = np.asarray(BinErrs)

    # write gaus parameters to file if selected
    if bWriteParsToFile:
        sPathToGaussSaveFile = sPathToFileSave.split('.')[0] + '.txt'
        fout = open(sPathToGaussSaveFile, 'w')
        for i in range(len(popts)):
            fout.write(str(BinCenters[i]) + '\t' + str(BinErrs[i]))
            for j in range(len(popts[0])):
                fout.write('\t' + str(popts[i][j]) + '\t' + str(pcovs[i][j][j]**0.5))
            fout.write('\n')
        fout.close()


    return (BinCenters, BinErrs, popts, pcovs)


def GetPreliminaryFit(x, y, xBins, qs=[16,50,84], FitType='Exp', BinsOmitUp=5, BinsOmitLow=5):
    BinCenters, Percentiles = GetBinPercentiles(x, y, xBins, qs=qs)

    AvgErr = 0.5 * (Percentiles.T[2] - Percentiles.T[0])

    if FitType == 'Exp':
        popt, pcov = curve_fit(FitExponential, BinCenters[BinsOmitUp:-BinsOmitLow],
                               Percentiles.T[1][BinsOmitUp:-BinsOmitLow], [1e5, 400],
                               sigma=AvgErr[BinsOmitUp:-BinsOmitLow])

    return popt


def GetLifetimeLogPlot(data, dataOmit, MaxValueDt, S2Min=5000, S2Max=200000):
    fig = plt.figure(figsize=(15, 8))
    Bins = [np.linspace(0,MaxValueDt,100), np.logspace(np.log10(S2Min), np.log10(S2Max), 100)]
#    xExcluded = data[Xcuts & (Xs2dt==0)]['dt']
#    yExcluded = data[Xcuts & (Xs2dt==0)]['s2']
#    xIncluded = data[Xcuts & Xs2dt]['dt']
#    yIncluded = data[Xcuts & Xs2dt]['s2']
    xExcluded = dataOmit.dt
    yExcluded = dataOmit.s2
    xIncluded = data.dt
    yIncluded = data.s2

    ExcludedHist, xBinEdges, yBinEdges = np.histogram2d(xExcluded, yExcluded, Bins)
    IncludedHist, xBinEdges, yBinEdges = np.histogram2d(xIncluded, yIncluded, Bins)
    X, Y = np.meshgrid(xBinEdges, yBinEdges)
    plt.pcolormesh(X, Y, ExcludedHist.T, norm=colors.LogNorm(), cmap='viridis', alpha=0.5)
    plt.pcolormesh(X, Y, IncludedHist.T, norm=colors.LogNorm(), cmap='viridis')

    return fig


def InsertRunsToDatasetList(RunNames, sPathToDatasetsList):
    FilesInList = ReadTotalRunsList(sPathToDatasetsList)
    for RunName in RunNames:
        if RunName in FilesInList:
            continue
        FilesInList.append(RunName)

    OutputFilelist(FilesInList, sPathToDatasetsList)
    return


def ReadElectronLifetimes(sPathToLifetimes):
    Unixtimes, UnixErrs, Lifetimes, LifetimeErrs = [], [], [], []
    try:
        fin = open(sPathToLifetimes, 'r')
        lines = fin.readlines()
        fin.close()
     
        for line in lines:
            contents = line[:-1].split('\t\t')
            if len(contents) != 4:
                continue
            Unixtimes.append(float(contents[0]))
            UnixErrs.append(float(contents[1]))
            Lifetimes.append(float(contents[2]))
            LifetimeErrs.append(float(contents[3]))

    except:
        print(sPathToLifetimes + ' does not exist')

    return (Unixtimes, UnixErrs, Lifetimes, LifetimeErrs)




def OutputElectronLifetimes(Unixtimes, UnixErrs, Lifetimes, LifetimeErrs, sPathToLifetimes):
    fout = open(sPathToLifetimes, 'w')

    for Unixtime, UnixErr, Lifetime, LifetimeErr in zip(Unixtimes, UnixErrs, Lifetimes, LifetimeErrs):
        fout.write(str(Unixtime) + '\t\t')
        fout.write(str(UnixErr) + '\t\t')
        fout.write(str(Lifetime) + '\t\t')
        fout.write(str(LifetimeErr) + '\n')
    fout.close()
    return


def GetInsertPosition(Unixtime, Unixtimes):
    if len(Unixtimes) == 0:
        return 0
    elif Unixtime < Unixtimes[0]:
        return 0
    elif Unixtime > Unixtimes[-1]:
        return len(Unixtimes)
    else:
        for i in range(len(Unixtimes)):
            if i == 0:
                continue
            if Unixtime > Unixtimes[i-1] and Unixtime < Unixtimes[i]:
                return i



def InsertLifetimeToFile(ValuesToInput, sPathToLifetimeFile, bOverwrite=False):
    Unixtime, UnixErr, Lifetime, LifetimeErr = ValuesToInput
    Unixtimes, UnixErrs, Lifetimes, LifetimeErrs = ReadElectronLifetimes(sPathToLifetimeFile)

    # if unixtime already in file, then see if want to overwrite
    if Unixtime in Unixtimes:
        if not bOverwrite:
            print('Unixtime already in ' + sPathToLifetimeFile.replace('/home/zgreene/', '') + ', returning without overwriting')
            return 
        else:
            print('Unixtime already in ' + sPathToLifetimeFile.replace('/home/zgreene/', '') + ', overwriting')
            Index = Unixtimes.index(Unixtime)
            Unixtimes[Index] = Unixtime
            UnixErrs[Index] = UnixErr
            Lifetimes[Index] = Lifetime
            LifetimeErrs[Index] = LifetimeErr

    # otherwise find where to insert such that chronological
    else:
        InsertPosition = GetInsertPosition(Unixtime, Unixtimes)
        Unixtimes.insert(InsertPosition, Unixtime)
        UnixErrs.insert(InsertPosition, UnixErr)
        Lifetimes.insert(InsertPosition, Lifetime)
        LifetimeErrs.insert(InsertPosition, LifetimeErr)

    OutputElectronLifetimes(Unixtimes, UnixErrs, Lifetimes, LifetimeErrs, sPathToLifetimeFile)
    return
