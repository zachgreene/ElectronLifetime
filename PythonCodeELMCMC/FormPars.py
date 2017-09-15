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
#    ChangeVal = -0.000144914181253
#    ChangeValErr = 1.20719462566e-05
#    ChangeVal = -0.000134390233
#    ChangeValErr = 0.000015843214
    # below values are generated from correcting PoRn alphas to Kr using prediction from early July
    ChangeVal = -0.000156713309
    ChangeValErr = 0.000015058217

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
        [1.90294571e+02, 1000.], # GXe volume outgassing, in unit of kg/day
        [1.94252374e+02, 1000.], # LXe volume outgassing, in unit of kg/day
        PurityDrops['unixtimes'], # time for the impurity change, after correction
            # 1496685600 from https://xenon-elog.lngs.infn.it/elog/XENON1T/484
        PurityDrops['types'],
        PurityDrops['values'],
        [[1471880000, 1472800000, -100.]],
#        1000., # GXe outgassing linear decreasing constant, in days.
#        1000., # LXe outgassing linear decreasing constant, in days.
        [[1471880000, 1472800000, 1.]], # fraction of GXe outgassing during period
        [[1475180000, 1475680000, 0.98]], # fraction of LXe outgassing during period
#        [1480144349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
        [1480344349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
        [1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#        [1500336000, 1501236000, 0.98], # period when getter is suspected to have lowered efficiency, July 18-28, 2017
        [1499731200, 1505230443, 0.98], # period when getter is suspected to have lowered efficiency, July 11-now, 2017
        ]

#    PurityDrops = GetPurityDrops()
#    default_pars = [
#        3.7e-3, # attaching rate from literature
#        7.12997918e+03, # initial GXe concentration
#        5.09257628e+01, # initial LXe concentration
#        4.02212420e-01, # impurity attaching prob for vaporization
#        4.10390827e-01, # impurity attaching prob for condensation
#        1.90294571e+02, # GXe volume outgassing, in unit of kg/day
#        1.94252374e+02, # LXe volume outgassing, in unit of kg/day
#        PurityDrops['unixtimes'], # time for the impurity change, after correction
#            # 1496685600 from https://xenon-elog.lngs.infn.it/elog/XENON1T/484
#        PurityDrops['types'],
#        PurityDrops['values'],
#        [[1471880000, 1472800000, -100.]],
#        1000., # GXe outgassing linear decreasing constant, in days.
#        1000., # LXe outgassing linear decreasing constant, in days.
#        [[1471880000, 1472800000, 1.]], # fraction of GXe outgassing during period
#        [[1475180000, 1475680000, 0.98]], # fraction of LXe outgassing during period
##        [1480144349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#        [1480344349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#        [1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
#        [1500336000, 1501236000, 0.98], # period when getter is suspected to have lowered efficiency, July 18-28, 2017
#        ]

    return default_pars


#def FormPars(x, MinUnixtime=-1, MaxUnixtime=-1):
def FormPars(x):
#    DaysSinceStart = 0
#    if MinUnixtime != -1:
#        DaysSinceStart = (MaxUnixtime - MinUnixtime) / 24. / 3600.
#        print('unixtime')

    PurityDrops = GetPurityDrops()
    pars = GetDefaultPars()
    NumPurityDrops = len(PurityDrops['unixtimes'])
    NumGXeOutgassingTerms = len(pars[5])
    NumLXeOutgassingTerms = len(pars[6])
    NumOutgassingTerms = NumGXeOutgassingTerms + NumLXeOutgassingTerms
    NumTerms = NumPurityDrops + NumOutgassingTerms
    
    if len(x)<13+NumPurityDrops:
#    	return GetDefaultPars()
        return pars
#    print("x=")
#    print(x)
    IfOutOfBoundary = False
#    print(NumTerms)
    for i, y in enumerate(x):
#    	if y<0 and (not i==11):
#    		IfOutOfBoundary = True
#    	if i==11 and y>0:
#    		IfOutOfBoundary = True
#    	if (i==2 or i==3 or i==15 or i==16 or i==17) and y>1:
#    		IfOutOfBoundary = True
#        if (i<4) and (y<0):
#            IfOutOfBoundary = True
#            print(1, i, y)
#        if (i==2 or i==3 or i>=7+NumTerms) and y>1:
#            IfOutOfBoundary = True
#            print(2, i, y)
#        elif (i==4) and (NumGXeOutgassingTerms > 0):
#        elif (i>=4) and (i<(4+NumGXeOutgassingTerms)):
#            if (i==4) and (y<0):
#            if (y<0):
#                IfOutOfBoundary = True
#            if i>4 and y<DaysSinceStart: 
#                IfOutOfBoundary = True
#                print('gxe out')
#                print(3, i, y)
#        elif (NumLXeOutgassingTerms > 0):
#        elif (i>=(4+NumGXeOutgassingTerms)) and (i<(4+NumGXeOutgassingTerms+NumLXeOutgassingTerms)):
#            if (i==(4+NumGXeOutgassingTerms)) and (y<0):
#            if (y<0):
#                IfOutOfBoundary = True
#            if (i>4+NumGXeOutgassingTerms) and y<DaysSinceStart: 
#                IfOutOfBoundary = True
#                print('lxe out')
#                print(4, i, y)
#        elif (i==(4+NumTerms)):
#            if y>0:
#                IfOutOfBoundary = True
#                print(5, i, y)
#        elif y<0:
#            IfOutOfBoundary = True
#            print(6, i, y)

    	if y<0 and (not i==4+NumTerms):
    		IfOutOfBoundary = True
    	if i==4+NumTerms and y>0:
    		IfOutOfBoundary = True
    	if (i==2 or i==3 or i>=7+NumTerms) and y>1:
    		IfOutOfBoundary = True
#    pars = GetDefaultPars()
    pars[1] = x[0] # initial GXe concentration
    pars[2] = x[1] # initial LXe concentration
    pars[3] = x[2] # vaporization attaching prob
    pars[4] = x[3] # condensation attaching prob
    for i in range(NumGXeOutgassingTerms):
        pars[5][i] = x[4+i]
    for i in range(NumLXeOutgassingTerms):
        pars[6][i] = x[4+i+NumLXeOutgassingTerms]
#    pars[5] = x[4] # GXe outgassing
#    pars[6] = x[5] # LXe outgassing
#    print(x[6:11])
    for i in range(NumPurityDrops):
        pars[9][i] = x[4+i+NumOutgassingTerms]
#        print(x[6+i], pars[9][i])
    pars[10][0][2] = x[4+NumTerms] # the amount of outgassing in GXe changing due to gas-only flow
#    pars[9]       = x[5+NumTerms] # GXe outgassing exponential decreasing constant, in days.
#    pars[10]      = x[6+NumTerms] # LXe outgassing exponential decreasing constant, in days.
    pars[11][0][2] = x[5+NumTerms] # fraction of GXe outgassing during gas-only circulation
    pars[12][0][2] = x[6+NumTerms] # fraction of LXe outgassing during PUR upgrade
    pars[13][2]    = x[7+NumTerms] # lowered efficiency
    pars[14][2]    = x[8+NumTerms] # lowered efficiency for Rn calibration during Christmas
    pars[15][2]    = x[9+NumTerms] # lowered efficiency for end of July 2017
    
    return (pars, IfOutOfBoundary)


def GetInitialParametersMCMC():
    x0 = np.array([
        5.82212635e+03, # initial GXe concentration
        6.55483847e+01, # initial LXe concentration
        7.04240787e-01, # impuritiy attachment prob for vaporization
        2.20849318e-01, # impuritiy attachment prob for condensation
        4.09610831e+02, # GXe volume outgassig, in units of kg/day
        1000.,          # GXe outgassing linear decreasing constant, in days
#        5000.,          # GXe outgassing linear decreasing constant, in days
        6.34534001e+01, # LXe volume outgassig, in units of kg/day
        1000.,          # LXe outgassing linear decreasing constant, in days
#        5000.,          # GXe outgassing linear decreasing constant, in days
        1.00e-5,        # impurity change from power event
        3.77e-6,        # impurity change from LN2 test July 15, 2016
        3.77e-6,        # impurity change from power glitch Nov. 21, 2016
        1.77e-6,        # impurity change after earthquake late January, 2017
        1.77e-6,        # impurity change from gate washing June 5, 2017
        -90,            # amount of GXe outgassing from gas-only flow
                        # 1000.,
                        # 1000.,
        1.,             # fraction of outgassing during gas-only circulation
        0.5,            # fraction of outgassing during PUR upgrade
        0.5,            # lowered efficiency Nov. 28 to Dec. 6
                        # 0.2,
        0.5,            # lowered efficiency during Rn calibration for SR0, end of Dec. 2016
        0.5,            # lowered efficiency for end of July 2017
		])
    x0_steps = np.array([
        2e3,  # initial GXe concentration
        30.,  # initial LXe concentration
        0.1,  # impuritiy attachment prob for vaporization
        0.06, # impuritiy attachment prob for condensation
        150., # GXe volume outgassig, in units of kg/day
        300., # GXe outgassing linear decreasing constant, in days
#        3000., # GXe outgassing linear decreasing constant, in days
        20,   # LXe volume outgassig, in units of kg/day
        300., # LXe outgassing linear decreasing constant, in days
#        3000., # GXe outgassing linear decreasing constant, in days
        3e-6, # impurity change from power event
        1e-6, # impurity change from LN2 test July 15, 2016
        1e-6, # impurity change from power glitch Nov. 21, 2016
        2e-7, # impurity change after earthquake late January, 2017
        2e-7, # impurity change from gate washing June 5, 2017
        30,   # amount of GXe outgassing from gas-only flow
              # 1000.,
              # 1000.,
        0.3,  # fraction of outgassing during gas-only circulation
        0.2,  # fraction of outgassing during PUR upgrade
        0.2,  # lowered efficiency Nov. 28 to Dec. 6
              # 0.2,
        0.2,  # lowered efficiency during Rn calibration for SR0, end of Dec. 2016
        0.2,  # lowered efficiency for end of July 2017
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
                ['$\\tau_{\Lambda_g}$', 'GXe outgassing linear constant', 'days'],
                ['$\\tau_{\Lambda_l}$', 'LXe outgassing linear constant', 'days'],
                ['$\\tau_{\Lambda_g}$', 'GXe outgassing quadratic constant', 'days'],
                ['$\\tau_{\Lambda_l}$', 'LXe outgassing quadratic constant', 'days'],
                ['$\Delta I_1$', 'Impurity change during power event, June 14', 'mol'],
                ['$\Delta I_2$', 'Impurity change during LN2 test, July 15', 'mol'],
                ['$\Delta I_3$', 'Impurity change during power glitch, Nov. 21', 'mol'],
                ['$\Delta I_4$', 'Impurity change after earthquake, late January', 'mol'],
                ['$\Delta I_5$', 'Impurity change from gate washing, 17/06/05', 'mol'],
                ['$\Delta \Lambda_g$', 'Decrease in GXe outgassing from GXe-only flow', 'kg*ppb/day'],
#                ['$\\tau_{\Lambda_g}$', 'GXe outgassing linear constant', 'days'],
#                ['$\\tau_{\Lambda_l}$', 'LXe outgassing linear constant', 'days'],
                ['$f_l$', 'fraction of GXe outgassing during GXe-only circulation', ''],
                ['$f_g$', 'fraction of LXe outgassing during PUR upgrade', ''],
                ['$\\alpha_1 $', 'lowered getter efficiency, Nov. 28 - Dec. 5 2016', ''],
                ['$\\alpha_2 $', 'lowered getter efficiency, Dec. 19 - Dec 21 2016', ''],
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
        ImpactfulUnixtimes.append(default_pars[11][0][0])
        ImpactfulUnixtimes.append(default_pars[11][0][1])
        ImpactfulUnixtimes.append(default_pars[12][0][0])
        ImpactfulUnixtimes.append(default_pars[12][0][1])
        ImpactfulUnixtimes.append(default_pars[13][0])
        ImpactfulUnixtimes.append(default_pars[13][1])
        ImpactfulUnixtimes.append(default_pars[14][0])
        ImpactfulUnixtimes.append(default_pars[14][1])
        ImpactfulUnixtimes.append(default_pars[15][0])
        ImpactfulUnixtimes.append(default_pars[15][1])

        return ImpactfulUnixtimes
