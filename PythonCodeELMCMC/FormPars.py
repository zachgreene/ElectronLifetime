import numpy as np
import datetime

def GetMinTimeStamp():
    return '05/17/16 00:00:00 '

def GetMaxTimeStamp(DaysAfterLastPoint=30):
#        return '08/21/17 00:00:00 '
        d = datetime.datetime.now() + datetime.timedelta(DaysAfterLastPoint)
        return d.strftime('%m/%d/%y %H:%M:%S ')

def GetS1ExponentialConstant():
    return 2040.6

# function to correct lifetime values found with pax_v6.2.0 and earlier
# from https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=greene:update_electron_lifetime_model
def GetLifetimeCorrectionPAX():
    CorrectBeforeTime = 1478000000.
    LifetimeCorrection = 6.72734759404e-05
    LifetimeCorrectionUncertainty = 1.56664061785e-06

    return (CorrectBeforeTime, LifetimeCorrection, LifetimeCorrectionUncertainty)


# if was calculated before this time, slightly off because of PoRn combined fit
def GetPoRnCorrection():
    CorrectBeforeTime = 1486124188.
    LifetimeCorrection = 2.19733154993e-05
    LifetimeCorrectionUncertainty = 7.67416766802e-06

    return (CorrectBeforeTime, LifetimeCorrection, LifetimeCorrectionUncertainty)


def GetKrCorrection():
    ChangeVal = -0.000144914181253
    ChangeValErr = 1.20719462566e-05

    return (ChangeVal, ChangeValErr)


def GetCathodeVoltages():
    CathodeVoltages = [
            [[0, 1473956519], [10., 20.]],
            [[1473956519, 1473997763], [15., 15.]],
            [[1473997763, 1475301507], [10., 20.]],
            [[1475301507, 1475391507], [12., 12.]],
            [[1475391507, 1484768041], [10., 20.]],
            [[1484768041, 1484942279], [0., 0.]],
            [[1484942279, 1485445141], [8., 9.]],
            [[1485445141, 1485802500], [0., 0.]],
            [[1485802500, 1486054320], [7., 7.]],
            [[1486054320, 1487265420], [8., 8.]],
            [[1487265420, int(2**32-1)], [7., 10.]]
            ]

    for i in range(len(CathodeVoltages)):
        if i == 0:
            continue
        assert CathodeVoltages[i][0][0] == CathodeVoltages[i-1][0][1]

    return CathodeVoltages

def GetSpecialPeriods():
    SpecialPeriods = [
                        [1465913920, 1466517600]
                        ]
    return SpecialPeriods

def GetGasOnlyPeriods():
    GasOnlyPeriods = [
                        [1471900000, 1472840000],
                        ]
    return GasOnlyPeriods

def GetMaximumHeatingPower():
    return 260.

def GetPurityDrops():
    PurityDrops = dict(
                    unixtimes = [1465937520,  # the amount impurity changed during power event
                                1468597800,   # the amount impurity changed during LN2 test @ July 15, 2016
                                1479772379,   # amount impurity changed during power glitch @ Nov. 21, 2016`
                                1485951100,   # the amount impurity after earthquake in late January 2017
                                1496685600],  # amount of impurity from gate washing on June 5, 2017
                    types = [0, 0, 0, 1, 1],
                    values = [1.00613537e-04, 3.36014333e-05, 4.0e-6, 1.0e-6, 1.e-6]
                    )
    assert len(PurityDrops['unixtimes']) == len(PurityDrops['types'])
    assert len(PurityDrops['unixtimes']) == len(PurityDrops['values'])

    return PurityDrops


def GetScienceRunUnixtimes():
    ScienceRunUnixtimes = dict(
                            SR0 = [1479772800, 1484731512],
                            SR1 = [1486054320, 2486054320]
                            )
#    ScienceRunUnixtimes = [
#                            [1479772800, 1484731512],
#                            [1486054320, 2486054320]
#                            ]
    return ScienceRunUnixtimes



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
# 170403: 14 parameters
# 170405: 16 parameters
def GetDefaultPars():
    PurityDrops = GetPurityDrops()
    default_pars = [
        3.7e-3, # attaching rate from literature
        7.12997918e+03, # initial GXe concentration
        5.09257628e+01, # initial LXe concentration
        4.02212420e-01, # impurity attaching prob for vaporization
        4.10390827e-01, # impurity attaching prob for condensation
        1.90294571e+02, # GXe volume outgassing, in unit of kg/day
        1.94252374e+02, # LXe volume outgassing, in unit of kg/day
        PurityDrops['unixtimes'], # time for the impurity change, after correction
            # 1496685600 from https://xenon-elog.lngs.infn.it/elog/XENON1T/484
        PurityDrops['types'],
        PurityDrops['values'],
        [[1471880000, 1472800000, -100.]],
        1000., # GXe outgassing linear decreasing constant, in days.
        1000., # LXe outgassing linear decreasing constant, in days.
        [[1471880000, 1472800000, 1.]], # fraction of GXe outgassing during period
        [[1475180000, 1475680000, 0.98]], # fraction of LXe outgassing during period
#        [1480144349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
        [1480344349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
        [1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
        [1500336000, 1501236000, 0.98], # period when getter is suspected to have lowered efficiency, July 18-28, 2017
        ]

    return default_pars


def FormPars(x):
    PurityDrops = GetPurityDrops()
    NumPurityDrops = len(PurityDrops['unixtimes'])
    
#    if len(x)<13+NumPurityDrops:
    if len(x)<13+NumPurityDrops:
    	return GetDefaultPars()
    print("x=")
    print(x)
    IfOutOfBoundary = False
    for i, y in enumerate(x):
    	if y<0 and (not i==11):
    		IfOutOfBoundary = True
    	if i==11 and y>0:
    		IfOutOfBoundary = True
    	if (i==2 or i==3 or i==15 or i==16 or i==17) and y>1:
    		IfOutOfBoundary = True
    pars = GetDefaultPars()
    pars[1] = x[0] # initial GXe concentration
    pars[2] = x[1] # initial LXe concentration
    pars[3] = x[2] # vaporization attaching prob
    pars[4] = x[3] # condensation attaching prob
    pars[5] = x[4] # GXe outgassing
    pars[6] = x[5] # LXe outgassing
#    print(x[6:11])
    
    for i in range(NumPurityDrops):
        pars[9][i] = x[6+i]
#        print(x[6+i], pars[9][i])
    pars[10][0][2] = x[6+NumPurityDrops] # the amount of outgassing in GXe changing due to gas-only flow
    pars[11] = x[7+NumPurityDrops] # GXe outgassing exponential decreasing constant, in days.
    pars[12] = x[8+NumPurityDrops] # LXe outgassing exponential decreasing constant, in days.
    pars[13][0][2] = x[9+NumPurityDrops] # fraction of GXe outgassing during gas-only circulation
    pars[14][0][2] = x[10+NumPurityDrops] # fraction of LXe outgassing during PUR upgrade
    pars[15][2] = x[11+NumPurityDrops] # lowered efficiency
    pars[16][2] = x[12+NumPurityDrops] # lowered efficiency for Rn calibration during Christmas
    pars[17][2] = x[13+NumPurityDrops] # lowered efficiency for end of July 2017
    
    return (pars, IfOutOfBoundary)


def GetInitialParametersMCMC():
    x0 = np.array([
        5.82212635e+03,
        6.55483847e+01,
        7.04240787e-01,
        2.20849318e-01,
        4.09610831e+02,
        6.34534001e+01,
        1.00e-5,
        3.77e-6,
        3.77e-6,
        1.77e-6,
        1.77e-6,
#		2.77e-6,
#		2.77e-6,
#		1.77e-6,
        -90,
        1000.,
        1000.,
        1.,
#		0.9,
        0.5,
#		0.95,
        0.5,
#		0.2,
        0.5,
        0.5,
		])
    x0_steps = np.array([
#		5e3,
#		3e3,
        2e3,
#		20.,
        30.,
        0.2,
#		0.3,
        0.06,
        150.,
        20,
        3e-6,
        1e-6,
        1e-6,
#		1e-6,
        2e-7,
        2e-7,
        30,
        300.,
        300.,
        0.3,
        0.2,
        0.2,
        0.2,
        0.2,
        ])

    return (x0, x0_steps)


def GetParInfo():
    ParInfo = [
                ['$I_g^0$', 'initial GXe concentration', 'ppb'],
                ['$I_l^0$', 'initial LXe concentration', 'ppb'],
                ['$\epsilon_1$', '$O_2$ attach probability for LXe vaporization', ''],
                ['$\epsilon_2$', '$O_2$ attach probability for GXe condensation', ''],
                ['$\Lambda_g^0$', 'initial GXe outgassing rate', 'kg*ppb/day'],
                ['$\Lambda_l^0$', 'initial LXe outgassing rate', 'kg*ppb/day'],
                ['$\Delta I_1$', 'Impurity change during power event, June 14', 'mol'],
                ['$\Delta I_2$', 'Impurity change during LN2 test, July 15', 'mol'],
                ['$\Delta I_3$', 'Impurity change during power glitch, Nov. 21', 'mol'],
                ['$\Delta I_4$', 'Impurity change after earthquake, late January', 'mol'],
                ['$\Delta I_5$', 'Impurity change from gate washing, 17/06/05', 'mol'],
                ['$\Delta \Lambda_g$', 'Decrease in GXe outgassing from GXe-only flow', 'kg*ppb/day'],
                ['$\\tau_{\Lambda_g}$', 'GXe outgassing linear constant', 'days'],
                ['$\\tau_{\Lambda_l}$', 'LXe outgassing linear constant', 'days'],
                ['$f_l$', 'fraction of GXe outgassing during GXe-only circulation', ''],
                ['$f_g$', 'fraction of LXe outgassing during PUR upgrade', ''],
                ['$\\alpha_1 $', 'lowered getter efficiency, Nov. 28 - Dec. 5 2016', ''],
                ['$\\alpha_2 $', 'lowered getter efficiency, Dec. 19 - Dec 21 2016', '']
                ['$\\alpha_3 $', 'lowered getter efficiency, July 18 - July 28 2017', '']
                ]

    return ParInfo


def GetImpactfulUnixtimes():
        ScienceRunUnixtimes = GetScienceRunUnixtimes()
        CathodeVoltages = GetCathodeVoltages()
        default_pars = GetDefaultPars()

        ImpactfulUnixtimes = []
        for ScienceRunUnixtimes in ScienceRunUnixtimes.values():
            for ScienceRunUnixtime in ScienceRunUnixtimes:
                if ScienceRunUnixtime > 2e9:
                    continue
                ImpactfulUnixtimes.append(ScienceRunUnixtime)

        for CathodeVoltage in CathodeVoltages:
            if CathodeVoltage[0][0] != 0:
                ImpactfulUnixtimes.append(CathodeVoltage[0][0])

        ImpactfulUnixtimes.append(default_pars[7][0])
        ImpactfulUnixtimes.append(default_pars[7][1])
        ImpactfulUnixtimes.append(default_pars[7][2])
        ImpactfulUnixtimes.append(default_pars[7][3])
        ImpactfulUnixtimes.append(default_pars[7][4])
        ImpactfulUnixtimes.append(default_pars[10][0][0])
        ImpactfulUnixtimes.append(default_pars[10][0][1])
        ImpactfulUnixtimes.append(default_pars[13][0][0])
        ImpactfulUnixtimes.append(default_pars[13][0][1])
        ImpactfulUnixtimes.append(default_pars[14][0][0])
        ImpactfulUnixtimes.append(default_pars[14][0][1])
        ImpactfulUnixtimes.append(default_pars[15][0])
        ImpactfulUnixtimes.append(default_pars[15][1])
        ImpactfulUnixtimes.append(default_pars[16][0])
        ImpactfulUnixtimes.append(default_pars[16][1])
        ImpactfulUnixtimes.append(default_pars[17][0])
        ImpactfulUnixtimes.append(default_pars[17][1])

        return ImpactfulUnixtimes
