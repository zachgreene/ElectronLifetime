import pandas as pd
import numpy as np

def GetInitialCuts(pax_version, RunSources):
    cuts = dict(
            LowLimitS1 = 30000,
            UpLimitS1 = 160000,
            LowLimitS2 = 20000,
            UpLimitS2 = 200000,
            LowLimitAsymS1 = 0.02,
            UpLimitAsymS1 = 0.5,
            LowLimitAsymS2 = 0.5,
            UpLimitAsymS2 = 0.75,
            RadialLimit = (1000.)**0.5,
            LowBoundZ = -100,
            UpBoundZ = 0,
            MaxValueDt = 730,
            NumBinsLifetime = 20,
            BinsOmitUp = 1,
            BinsOmitLow = 1,
            FactorBelow = 0.6,
            FactorAbove = 1.5
           )

    # change limits for earlier pax version
    if pax_version == '6.2.0':
        cuts['LowLimitS1'] = 20000
        cuts['UpLimitS1'] = 70000
        cuts['LowLimitS2'] = 15000
        cuts['UpLimitS2'] = 125000

    if int(pax_version.replace('.', '')) >= 665:
        cuts['LowLimitS1'] = 25000
        cuts['UpLimitS1'] = 80000
        cuts['LowLimitS2'] = 18000
        cuts['UpLimitS2'] = 200000
        cuts['BinsOmitLow'] = 2

    # rn220 has large amount of low-z events
    if 'Rn220' in RunSources:
        cuts['BinsOmitLow'] = 2
        cuts['NumBinsLifetime'] = 40,

    return cuts


def SetS2DtLow(data, PrelimOpts, FactorBelow):
    p0, p1 = PrelimOpts
    Xeval = (data['s2'] > FactorBelow*p0*np.exp(-data['dt']/p1))
    return Xeval


def SetS2DtUp(data, PrelimOpts, FactorAbove):
    p0, p1 = PrelimOpts
    Xeval = (data['s2'] < FactorAbove*p0*np.exp(-data['dt']/p1))
    return Xeval

def EvalCut(data, cuts, CutToEval):
    Xeval = []
    if CutToEval == 'Xs1':
        Xeval = (data['s1'] > cuts['LowLimitS1']) & (data['s1'] < cuts['UpLimitS1'])
    elif CutToEval == 'Xs2':
        Xeval = (data['s2'] > cuts['LowLimitS2']) & (data['s2'] < cuts['UpLimitS2'])
    elif CutToEval == 'Xs1asym':
        Xeval = (data['s1_aft'] > cuts['LowLimitAsymS1']) & (data['s1_aft'] < cuts['UpLimitAsymS1'])
    elif CutToEval == 'Xs2asym':
        Xeval = (data['s2_aft'] > cuts['LowLimitAsymS2']) & (data['s2_aft'] < cuts['UpLimitAsymS2'])
    elif CutToEval == 'Xr':
        Xeval = (data['r'] < cuts['RadialLimit'])
    elif CutToEval == 'Xz':
        Xeval = (data['z'] < cuts['UpBoundZ']) & (data['z'] > cuts['LowBoundZ'])
    return Xeval
