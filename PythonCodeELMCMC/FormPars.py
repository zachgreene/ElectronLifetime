import numpy as np

def GetMinTimeStamp():
	return '05/17/16 00:00:00 '

def GetMaxTimeStamp():
	return '08/21/17 00:00:00 '

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
	default_pars = [
		3.7e-3, # attaching rate from literature
		7.12997918e+03, # initial GXe concentration
		5.09257628e+01, # initial LXe concentration
		4.02212420e-01, # impurity attaching prob for vaporization
		4.10390827e-01, # impurity attaching prob for condensation
		1.90294571e+02, # GXe volume outgassing, in unit of kg/day
		1.94252374e+02, # LXe volume outgassing, in unit of kg/day
		[1465937520, 1468597800, 1479772379, 1485951100, 1496685600], # time for the impurity change, after correction
                # 1496685600 from https://xenon-elog.lngs.infn.it/elog/XENON1T/484
		[0, 0, 0, 1, 1],
		[1.00613537e-04, 3.36014333e-05, 4.0e-6, 1.0e-6, 1.e-6],
		[[1471880000, 1472800000, -100.]],
		1000., # GXe outgassing linear decreasing constant, in days.
		1000., # LXe outgassing linear decreasing constant, in days.
		[[1471880000, 1472800000, 1.]], # fraction of GXe outgassing during period
		[[1475180000, 1475680000, 0.98]], # fraction of LXe outgassing during period
#		[1480144349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
		[1480344349, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
		[1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
		]

	return default_pars


def FormPars(x):
	if len(x)<18:
		return GetDefaultPars()
#	print("x=")
#	print(x)
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
	pars[9][0] = x[6] # the amount impurity changed during power event
	pars[9][1] = x[7] # the amount impurity changed during LN2 test @ July 15, 2016
	pars[9][2] = x[8] # the amount impurity changed during power glitch @ Nov. 21, 2016
	pars[9][3] = x[9] # the amount impurity after earthquake in late January 2017
	pars[9][4] = x[10] # amount of impurity from gate washing on June 5, 2017
	pars[10][0][2] = x[11] # the amount of outgassing in GXe changing due to gas-only flow
	pars[11] = x[12] # GXe outgassing exponential decreasing constant, in days.
	pars[12] = x[13] # LXe outgassing exponential decreasing constant, in days.
	pars[13][0][2] = x[14] # fraction of GXe outgassing during gas-only circulation
	pars[14][0][2] = x[15] # fraction of LXe outgassing during PUR upgrade
	pars[15][2] = x[16] # lowered efficiency
	pars[16][2] = x[17] # lowered efficiency for Rn calibration during Christmas
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
                ['$\\alpha_1 $', 'lowered getter efficiency, Nov. 28 - Dec. 5', ''],
                ['$\\alpha_2 $', 'lowered getter efficiency, Dec. 19 - Dec 21', '']
                ]

    return ParInfo


def GetScienceRunUnixtimes():
    ScienceRunUnixtimes = [
                            [1479772800, 1484731512]
                            ]
    return ScienceRunUnixtimes


def GetImpactfulUnixtimes():
        ScienceRunUnixtimes = GetScienceRunUnixtimes()
        default_pars = GetDefaultPars()

        ImpactfulUnixtimes = []
        for ScienceRun in ScienceRunUnixtimes:
            for StartEndTime in ScienceRun:
                ImpactfulUnixtimes.append(StartEndTime)

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

        return ImpactfulUnixtimes
