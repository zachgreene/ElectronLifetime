import re
import numpy as np
import pandas as pd
import scipy as scp
import os, sys
import time, datetime

import ROOT
from ROOT import TH1
from ROOT import TH1D
from ROOT import TGraph

import FormPars

def GetFrame(hname, MinUnixTime, MaxUnixTime, MinValue, MaxValue, xtitle, ytitle):
    hframe = TH1D(hname, "", 100, MinUnixTime/3600./24., MaxUnixTime/3600./24.)
    hframe.SetStats(0)
    hframe.GetXaxis().SetLabelFont(132)
    hframe.GetXaxis().SetTitleFont(132)
    hframe.GetYaxis().SetLabelFont(132)
    hframe.GetYaxis().SetTitleFont(132)
    hframe.GetXaxis().SetLabelSize(0.)
    hframe.GetXaxis().SetTitleSize(0.)
    hframe.GetYaxis().SetLabelSize(0.)
    hframe.GetYaxis().SetTitleSize(0.)
    hframe.GetYaxis().SetRangeUser(MinValue, MaxValue)
    if len(xtitle)>0:
        hframe.GetXaxis().SetLabelSize(0.06)
        hframe.GetXaxis().SetTitleSize(0.06)
        hframe.GetXaxis().CenterTitle()
        hframe.GetXaxis().SetTitle(xtitle)
    if len(ytitle)>0:
        hframe.GetYaxis().SetLabelSize(0.06)
        hframe.GetYaxis().SetTitleSize(0.06)
        hframe.GetYaxis().CenterTitle()
        hframe.GetYaxis().SetTitle(ytitle)
    return hframe

def GetFrameSecondUnit(hname, MinUnixTime, MaxUnixTime, MinValue, MaxValue, xtitle, ytitle):
    hframe = TH1D(hname, "", 100, MinUnixTime, MaxUnixTime)
    hframe.SetStats(0)
    hframe.GetXaxis().SetLabelFont(132)
    hframe.GetXaxis().SetTitleFont(132)
    hframe.GetYaxis().SetLabelFont(132)
    hframe.GetYaxis().SetTitleFont(132)
    hframe.GetXaxis().SetLabelSize(0.)
    hframe.GetXaxis().SetTitleSize(0.)
    hframe.GetYaxis().SetLabelSize(0.)
    hframe.GetYaxis().SetTitleSize(0.)
    hframe.GetYaxis().SetRangeUser(MinValue, MaxValue)
    if len(xtitle)>0:
        hframe.GetXaxis().SetLabelSize(0.06)
        hframe.GetXaxis().SetTitleSize(0.06)
        hframe.GetXaxis().CenterTitle()
        hframe.GetXaxis().SetTitle(xtitle)
    if len(ytitle)>0:
        hframe.GetYaxis().SetLabelSize(0.06)
        hframe.GetYaxis().SetTitleSize(0.06)
        hframe.GetYaxis().CenterTitle()
        hframe.GetYaxis().SetTitle(ytitle)
    return hframe


# Get unixtime from time stamp
def GetUnixTimeFromTimeStampTool(timestamp):
    if len(timestamp)<9:
        return None
    Contents = timestamp[0:9].split("_")
    #print(Contents)
    Year = "2016"
    Month = Contents[0][0:2]
    Day = Contents[0][2:4]
    #print(Day)
    Hour = Contents[1][0:2]
    Minute = Contents[1][2:4]
    Second = "0"
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second))
    return time.mktime(dt.timetuple())

# Get unixtime from time stamp
# but also return the datetime object
def GetUnixTimeAndDatetimeFromTimeStampTool(timestamp):
    if len(timestamp)<11:
        return (None, None)
    Contents = timestamp[0:11].split("_")
    #print(Contents)
    Year = "20"+Contents[0][0:2]
    Month = Contents[0][2:4]
    Day = Contents[0][4:6]
    #print(Day)
    Hour = Contents[1][0:2]
    Minute = Contents[1][2:4]
    Second = "0"
    #print(Year)
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second))
    return (dt,time.mktime(dt.timetuple()))

def GetUnixTimeFromTimeStamp(timestamp):
    #print(timestamp)
    reduce_stamp = timestamp[0:17]
    if len(timestamp)<18:
        return None
    Contents = reduce_stamp.split(" ")
    Month, Day, Year = Contents[0].split("/")
    Hour, Minute, Second = Contents[1].split(":")
    Year = "20"+Year
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second))
    return time.mktime(dt.timetuple())



# Get the dictionary containing 
# 1) unixtime
# 2) values
# from pandas series
#Translate one time stamp to unixtime
#yyyy-mm-dd HH:MM:SS.subsec
def TranslateTimestampToUnixtime(timestamp):
    TimeStamp = timestamp
    #print(TimeStamp)
    TimeStamp = TimeStamp.split(".")[0]
    Date, Hours = TimeStamp.split(" ")
    Year, Month, Day = Date.split("-")
    Hour, Minute, Second = Hours.split(":")
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second))
    return time.mktime(dt.timetuple())

#get the dictionary
# the series shall be pre-sorted already
def GetDictFromSeries(series):
    OutputData = {}
    Unixtimes = []
    Values = []
    for key, value in zip(series.keys(), series):
        # print(key)
        unixtime = TranslateTimestampToUnixtime(str(key))
        Unixtimes.append(unixtime)
        Values.append(value)
    OutputData['unixtimes'] = Unixtimes
    OutputData['values'] = Values
    return OutputData

# get region
def GetRegion(UnixTimes, LowerBoundaries, UpperBoundaries):
    if not len(LowerBoundaries)==len(UpperBoundaries):
        raise ValueError("Lists' length don't match")
    graph = TGraph()
    for i, (unixtime, lower) in enumerate(zip(UnixTimes,LowerBoundaries)):
        graph.SetPoint(i, unixtime, lower)
    for i, (unixtime, upper) in enumerate(zip(
                                                                      reversed(UnixTimes),
                                                                      reversed(UpperBoundaries)
                                                                     )):
        graph.SetPoint(i+len(UnixTimes), unixtime, upper)
    return graph

def PickupMCMCPosteriors(Samples, NumOfTrials):
    indexs = np.random.choice(len(Samples), NumOfTrials)
    OutputList = []
    for index in indexs:
        OutputList.append(Samples[index])
    return OutputList

# For Burn-In region
# count list depth
def depth(a):
    if type(a) is list or type(a) is np.ndarray:
        return 1+depth(a[0])
    return 0

# re-arrange the chain
# so that it become
# [Iteration][Walker_id][Par_id]
def ReshapeChain(chain, ndim, nwalkers, niterations):
    if not depth(chain)==3:
        raise ValueError("chain depth not standard!")
    OutputSamples = []
    # initialization
    # parameter loop
    for i in range(niterations):
        ListForAppend = []
        #walker loop
        for j in range(nwalkers):
            ListForAppend.append([])
        OutputSamples.append(ListForAppend)
    # fill the values
    # walker loop
    for i, SamplesOneWalker in enumerate(chain):
        #print(SamplesOneWalker)
        # iteration loop
        for k, SamplesOneIteration in enumerate(SamplesOneWalker):
            # Parameter loop
            for j, SampleOnePar in enumerate(SamplesOneIteration):
                #print(str(j)+"\t"+str(i)+"\t"+str(SampleOnePar))
                OutputSamples[k][i].append(SampleOnePar)
    return OutputSamples


# chain shall be a 3x3 vector
# [Par_id][Walker_id][Iteration]
def GetSamplesForPlot(chain, ndim, nwalkers, niterations):
    if not depth(chain)==3:
        raise ValueError("chain depth not standard!")
    OutputSamples = []
    # initialization
    # parameter loop
    for i in range(ndim):
        ListForAppend = []
        #walker loop
        for j in range(nwalkers):
            ListForAppend.append([])
        OutputSamples.append(ListForAppend)
    # fill the values
    # walker loop
    for i, SamplesOneWalker in enumerate(chain):
        #print(SamplesOneWalker)
        # iteration loop
        for SamplesOneIteration in SamplesOneWalker:
            # Parameter loop
            for j, SampleOnePar in enumerate(SamplesOneIteration):
                #print(str(j)+"\t"+str(i)+"\t"+str(SampleOnePar))
                OutputSamples[j][i].append(SampleOnePar)
    return OutputSamples

def GetBurnInCutoffSamples(chain, cutoff):
    OutputList = []
    #walker loop
    for SamplesOneWalker in chain:
        #print(SamplesOneWalker)
        # iteration loop
        for i, SamplesOneIteration in enumerate(SamplesOneWalker):
            if i<cutoff:
                continue
            OutputList.append(SamplesOneIteration)
    return OutputList

def GetBurnInCutoffSamplesV2(chain, cutofflower, cutoffupper):
    OutputList = []
    #walker loop
    for SamplesOneWalker in chain:
        #print(SamplesOneWalker)
        # iteration loop
        for i, SamplesOneIteration in enumerate(SamplesOneWalker):
            if i<cutofflower:
                continue
            if i>cutoffupper:
                continue
            OutputList.append(SamplesOneIteration)
    return OutputList


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


# function to correct lifetime values found with pax_v6.2.0 and earlier
# from https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=greene:update_electron_lifetime_model
def CorrectForPaxVersion(Lifetime, LifetimeErr):
    # found by plotting difference in inverse lifetimes
    LifetimeCorrectionFactor = 6.72734759404e-05
#    LifetimeCorrectionFactorUncertainty = 9.46946770525e-07
    LifetimeCorrectionFactorUncertainty = 1.56664061785e-06

    #idea is 1/tau_v6.2.0 - 1/tau_v6.4.0 = const_eff
    TrueLifetime = Lifetime/(1 - LifetimeCorrectionFactor*Lifetime)
    TrueLifetimeErr = LifetimeErr/(1. - LifetimeCorrectionFactor*Lifetime) + Lifetime/(1 - LifetimeCorrectionFactor*Lifetime)**2*(LifetimeCorrectionFactorUncertainty*Lifetime + LifetimeCorrectionFactor*LifetimeErr)

    return (TrueLifetime, TrueLifetimeErr)



def GetParameterPercentiles(samples, quantiles):
    ParameterQuantiles = np.percentile(samples, q=quantiles, axis=0)
    return ParameterQuantiles.T


def FormatNumberForWiki(aNumbers, bExp=False, PowersOfTen=None):
    aNumbers = np.asarray(aNumbers)
    if bExp:
        PowersOfTen = int(np.log10(np.abs(aNumbers[0])))
        if PowersOfTen < 0:
            PowersOfTen -= 1
        aNumbers /= 10**PowersOfTen

    sNumbers = ['%.2f' %(fNumber) for fNumber in aNumbers]

    StringForWiki = '%s_{-%s}^{+%s}' %(sNumbers[0], sNumbers[1], sNumbers[2])
    if bExp:
        StringForWiki += '\ e%i' %(PowersOfTen)

    return StringForWiki


def PrintParameterQuantilesWiki(samples, quantiles, ParNames, ParDescs, ParUnits):
    ParameterQuantiles = GetParameterPercentiles(samples, quantiles)
    ParErrorLow = ParameterQuantiles[:, 1] - ParameterQuantiles[:, 0]
    ParErrorUp = ParameterQuantiles[:, 2] - ParameterQuantiles[:, 1]

    print('^  Parameter  ^  Description  ^  Units  ^  Best Fit Value  ^')

    for i in range(samples.shape[1]):
#        if 'Delta' in ParNames[i]:
#            sParsToPrint = FormatNumberForWiki([ParameterQuantiles[i][1], ParErrorLow[i], ParErrorUp[i]], bExp=True)
#        else:
#            sParsToPrint = FormatNumberForWiki([ParameterQuantiles[i][1], ParErrorLow[i], ParErrorUp[i]])
        if 'Delta I' in ParNames[i]:
            StringForWiki = FormatNumberForWiki([ParameterQuantiles[i][1], ParErrorLow[i], ParErrorUp[i]], bExp=True)
        else:
            StringForWiki = FormatNumberForWiki([ParameterQuantiles[i][1], ParErrorLow[i], ParErrorUp[i]])
        print('|  ' + ParNames[i] + '  |  ', end='')
        print(ParDescs[i] + '  |  ', end='')
        print(ParUnits[i] + '  |  $', end='')
#        print(sParsToPrint[0] + '_{-', end='')
#        print(sParsToPrint[1] + '}^{+', end='')
        print(StringForWiki + '$  |')
#        print(str(ParameterQuantiles[i][1]) + '_{-', end='')
#        print(str(ParErrorLow[i]) + '}^{+', end='')
#        print(str(ParErrorUp[i]) + '}  |')

####################################
# Get Kr83m elifes 
####################################
def LoadLifetimesKr83(PathToFile='/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithKr83.txt'):
    fin = open(PathToFile)
    lines = fin.readlines()
    fin.close()
    
    KrUnixtimes = []
    KrUnixtimeErrors = []
    KrELifeValues = []
    KrELifeValueErrors = []
    
    for i, line in enumerate(lines):
        if line[0] == '#':
            continue
        contents = line[:-1].split(" ")
        unixtime = float(contents[0])
        value = float(contents[1])
        value_err = float(contents[2])
        KrUnixtimes.append(unixtime)
        KrUnixtimeErrors.append(0)
        KrELifeValues.append(value)
        KrELifeValueErrors.append(value_err)

    return (KrUnixtimes, KrUnixtimeErrors, KrELifeValues, KrELifeValueErrors)


def LoadLifetimesNorm(FitType, PathToFile):
    fin = open(PathToFile)
    lines = fin.readlines()
    fin.close()
    
    Unixtimes = []
    UnixtimeErrors = []
    ELifeValues = []
    ELifeValueErrors = []
    
    for i, line in enumerate(lines):
        if line[0] == '#':
            continue

        contents = line[:-1].split("\t\t")
        unixtime = float(contents[0])
        unixtime_err = float(contents[1])
        value = float(contents[2])
        value_err = float(contents[3])
        # account for LCE
        if FitType == 'SingleScatter':
            S1ExponentialConstant = FormPars.GetS1ExponentialConstant()
            value = value*S1ExponentialConstant / (S1ExponentialConstant - value)
            value_err = value * np.sqrt(np.power(value_err/value,2.0) +
                                        np.power(value_err/ (S1ExponentialConstant - value), 2.))
            if value_err < 0.:
                continue
        #if was calculated using pax_v6.2.0 or younger
        elif FitType == 'Rn' and unixtime < FormPars.GetLifetimeCorrectionPAX():
#            print(unixtime, value, value_err, *CorrectForPaxVersion(value, value_err))
            value, value_err = CorrectForPaxVersion(value, value_err)
        Unixtimes.append(unixtime)
        UnixtimeErrors.append(unixtime_err)
        ELifeValues.append(value)
        ELifeValueErrors.append(value_err)

    return (Unixtimes, UnixtimeErrors, ELifeValues, ELifeValueErrors)


def LoadFitData(FitType, PathToFile):
    if FitType == 'Kr83':
        return LoadLifetimesKr83()
    else:
        return LoadLifetimesNorm(FitType, PathToFile)
    return -1



def LoadPredictions(PredictionFile, LastPointUnixtime, DaysAfterLastPoint=0):
    fin = open(PredictionFile)
    lines = fin.readlines()
    fin.close()

    PredicionUnixtimes = []
    PredictedELifes = []
    PredictedELifeLows = []
    PredictedELifeUps = []
    PredictedELifeLowErrs = []
    PredictedELifeUpErrs = []
    for i, line in enumerate(lines):
        contents = line[:-1].split("\t\t")
        unixtime = float(contents[0])
        elife = float(contents[1])
        elife_lower = float(contents[2])
        elife_upper = float(contents[3])
        if unixtime > LastPointUnixtime + DaysAfterLastPoint*3600.*24.:
            break
        PredicionUnixtimes.append(unixtime)
        PredictedELifes.append(elife)
        PredictedELifeLows.append(elife_lower)
        PredictedELifeUps.append(elife_upper)
        PredictedELifeLowErrs.append( elife_lower - elife)
        PredictedELifeUpErrs.append(elife_upper - elife)
#        PredictedELifeLowErrs.append( (elife_lower - elife)/elife)
#        PredictedELifeUpErrs.append((elife_upper - elife)/elife)

    return (PredicionUnixtimes,
            PredictedELifes,
            PredictedELifeLows,
            PredictedELifeUps,
            PredictedELifeLowErrs,
            PredictedELifeUpErrs)


def ChangeElectronLifetime(ELifes, ELifeUps, ELifeLows, ELifeUpErrs, ELifeLowErrs, ChangeVal, ChangeValErr):

    NewLifetimes = []
    NewLifetimeUps = []
    NewLifetimeLows = []
    NewLifetimeUpErrs = []
    NewLifetimeLowErrs = []

    for ELife, ELifeUp, ELifeLow, ELifeUpErr, ELifeLowErr in zip(
                                        ELifes, ELifeUps, ELifeLows, ELifeUpErrs, ELifeLowErrs
                                        ):

        InverseNewLifetime = 1./ELife + ChangeVal  # as defined in GetRelationshipBetweenRnKr.py
        InverseNewLifetimeUp = 1./ELifeUp + ChangeVal
        InverseNewLifetimeLow = 1./ELifeLow + ChangeVal
        InverseNewLifetimeUpErr = -((ELifeUpErr/(ELifeUp**2))**2 + ChangeValErr**2)**0.5
        InverseNewLifetimeLowErr = ((ELifeLowErr/(ELifeLow**2))**2 + ChangeValErr**2)**0.5

        NewLifetimes.append(1./InverseNewLifetime)
        NewLifetimeUps.append(1./InverseNewLifetimeUp)
        NewLifetimeLows.append(1./InverseNewLifetimeLow)
        NewLifetimeUpErrs.append(InverseNewLifetimeUpErr/(InverseNewLifetimeUp**2.))
        NewLifetimeLowErrs.append(InverseNewLifetimeLowErr/(InverseNewLifetimeLow**2.))


    return (NewLifetimes,
            NewLifetimeUps,
            NewLifetimeLows,
            NewLifetimeUpErrs,
            NewLifetimeLowErrs)
