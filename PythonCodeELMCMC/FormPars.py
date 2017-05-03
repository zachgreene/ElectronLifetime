import numpy as np

def GetMinTimeStamp():
	return '05/17/16 00:00:00 '

def GetMaxTimeStamp():
	return '03/31/17 00:00:00 '

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

def GetDefaultPars():
	default_pars = [
		3.7e-3, # attaching rate from literature
		7.12997918e+03, # initial GXe concentration
		5.09257628e+01, # initial LXe concentration
		4.02212420e-01, # impurity attaching prob for vaporization
		4.10390827e-01, # impurity attaching prob for condensation
		1.90294571e+02, # GXe volume outgassing, in unit of kg/day
		1.94252374e+02, # LXe volume outgassing, in unit of kg/day
		[1465937520, 1468597800, 1479772379], # time for the impurity change, after correction
		[0, 0, 0],
		[1.00613537e-04, 3.36014333e-05, 1.0e-6],
		[[1471880000, 1472800000, -100.]],
		1000., # GXe outgassing linear decreasing constant, in days.
		1000., # LXe outgassing linear decreasing constant, in days.
		[[1487670000, 1000.], ], # additional LXe outgassing linear decreasing 
		[1480317149, 1480926700, 0.98], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
		[1482175745 - 2.*3600., 1482351960 + 2.*3600., 0.2], # periods when getter is suspected to have lowered efficiency, roughly from 11-28 to 12-06
		[1485351200, 1486272625, 0.9], # periods when getter is suspected to have lowered efficiency after earthquake
		]

	return default_pars


def FormPars(x):
	if len(x)<16:
		return GetDefaultPars()
	print("x=")
	print(x)
	IfOutOfBoundary = False
	for i, y in enumerate(x):
		if y<0 and (not i==9):
			IfOutOfBoundary = True
		if i==9 and y>0:
			IfOutOfBoundary = True
		if (i==2 or i==3 or i==13 or i==14 or i==15) and y>1:
			IfOutOfBoundary = True
	pars = GetDefaultPars()
	pars[1] = x[0] # initial GXe concentration
	pars[2] = x[1] # initial LXe concentration
	pars[3] = x[2] # vaporization attaching prob
	pars[4] = x[3] # condensation attaching prob
	pars[5] = x[4] # GXe outgassing
	pars[6] = x[5] # LXe outgassing
	pars[9][0] = x[6] # the amount impurity changed during power event
	pars[9][1] = x[7] # the amount impurity changed during LN2 test @ July 15th
	pars[9][2] = x[8] # the amount impurity changed during power glitch @ Nov. 21th
	pars[10][0][2] = x[9] # the amount of outgassing in GXe changing due to gas-only flow
	pars[11] = x[10] # GXe outgassing exponential decreasing constant, in days.
	pars[12] = x[11] # LXe outgassing exponential decreasing constant, in days.
	pars[13][0][1] = x[12] # LXe outgassing linear decreasing constant
	pars[14][2] = x[13] # lowered efficiency
	pars[15][2] = x[14] # lowered efficiency for Rn calibration during Christmas
	pars[16][2] = x[15] # lowered efficiency for after earthquake in January
	return (pars, IfOutOfBoundary)


def GetIniitalParametersMCMC():
	x0 = np.array([5.82212635e+03,
		6.55483847e+01,
		7.04240787e-01,
		2.20849318e-01,
		4.09610831e+02,
		6.34534001e+01,
		1.00e-5,
		2.77e-6,
		2.77e-7,
		-90,
		1000.,
		1000.,
		1000.,
		0.95,
		0.2,
		0.9,
		])
	x0_steps = np.array([
		5e2,
		6.,
		0.07,
		0.02,
		4.,
		0.6,
		1.e-6,
		2.77e-7,
		2.77e-8,
		4.5,
		1.,
		1.,
		1.,
		0.05,
		0.02,
		0.05,
		])

	return (x0, x0_steps)
