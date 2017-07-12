import numpy as np
import scipy as sp

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt
import time
import pickle
import sys
from Tools import *

if len(sys.argv)<2:
    print("======== Syntax: =======")
    print("python DrawPlotForElectronLifetimePrediction.py .....")
    print("< elife data txt file> ")
    print("<Rn elife data txt file>")
    print("<Kr83 elife data txt file>")
    print("< prediction txt file> ")
    print("< days to show after last data point >")
    print("<save fig name (rel.)>")
    exit()

ELifeDataFile = sys.argv[1]
RnELifeDataFile = sys.argv[2]
Kr83ELifeDataFile = sys.argv[3]
PredictionFile = sys.argv[4]
DaysAfterLastPoint = float(sys.argv[5])
FigureSaveName = sys.argv[6]

S1ExponentialConstant = 2040. # us. 
FitterUsedS1ExponentialConstant = 2040.

#######################################
### Get the elife data
#######################################
# electron lifetime
fin1 = open(ELifeDataFile)
lines1 = fin1.readlines()
fin1.close()

UnixTimes = []
UnixTimeErrors = []
ELifeValues = []
ELifeValueErrors = []
for i, line in enumerate(lines1):
    contents = line[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    # correct the e-life by the S1 term
    value = value*S1ExponentialConstant / (S1ExponentialConstant - value)
    value_err = value*np.sqrt( np.power(value_err/value,2.0)+np.power(value_err/ (S1ExponentialConstant - value), 2.) )
    UnixTimes.append(unixtime)
    UnixTimeErrors.append(unixtime_err)
    ELifeValues.append(value)
    ELifeValueErrors.append(value_err)
FirstPointUnixTime = UnixTimes[0]
LastPointUnixTime = UnixTimes[len(UnixTimes)-1]

######################################
## Get Rn elife data
######################################
fin3 = open(RnELifeDataFile)
lines3 = fin3.readlines()
fin3.close()
#print(lines3)

RnUnixtimes = []
RnUnixtimeErrors = []
RnELifeValues = []
RnELifeValueErrors = []

for i, line3 in enumerate(lines3):
    if line3[0] == '#':
#        print(line3)
        continue
    contents = line3[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    #if was calculated using pax_v6.2.0 or younger
#    if unixtime < 1478000000 or (unixtime > 1484900000 and unixtime < 1486100000):
    if unixtime < 1478000000:
        print(unixtime, value, value_err, *CorrectForPaxVersion(value, value_err))
        value, value_err = CorrectForPaxVersion(value, value_err)
    RnUnixtimes.append(unixtime)
    RnUnixtimeErrors.append(unixtime_err)
    RnELifeValues.append(value)
    RnELifeValueErrors.append(value_err)

LastPointUnixTime = RnUnixtimes[-1]

CutID = 0
for i, unixtime in enumerate(UnixTimes):
    if unixtime>RnUnixtimes[0]:
        CutID = i
        break
UnixTimes = UnixTimes[:CutID]
UnixTimeErrors = UnixTimeErrors[:CutID]
ELifeValues = ELifeValues[:CutID]
ELifeValueErrors = ELifeValueErrors[:CutID]
####################################
# Get Kr83m elifes 
####################################
fin4 = open(Kr83ELifeDataFile)
lines4 = fin4.readlines()
fin4.close()
#print(lines4)

KrUnixtimes = []
KrUnixtimeErrors = []
KrELifeValues = []
KrELifeValueErrors = []

for i, line4 in enumerate(lines4):
    contents = line4[:-1].split(" ")
    unixtime = float(contents[0])
    value = float(contents[1])
    value_err = float(contents[2])
    KrUnixtimes.append(unixtime)
    KrUnixtimeErrors.append(0)
    KrELifeValues.append(value)
    KrELifeValueErrors.append(value_err)


####################################
# Get Xe129m elifes 
####################################


Xe129mELifeDataFile = '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithXe129m.txt'
Xe131mELifeDataFile = '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithXe131m.txt'
fin5 = open(Xe129mELifeDataFile)
lines5 = fin5.readlines()
fin5.close()

Xe129mUnixtimes = []
Xe129mUnixtimeErrors = []
Xe129mELifeValues = []
Xe129mELifeValueErrors = []

for i, line5 in enumerate(lines5):
    contents = line5[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    Xe129mUnixtimes.append(unixtime)
    Xe129mUnixtimeErrors.append(unixtime_err)
    Xe129mELifeValues.append(value)
    Xe129mELifeValueErrors.append(value_err)


####################################
# Get Xe131m elifes 
####################################
fin6 = open(Xe131mELifeDataFile)
lines6 = fin6.readlines()
fin6.close()

Xe131mUnixtimes = []
Xe131mUnixtimeErrors = []
Xe131mELifeValues = []
Xe131mELifeValueErrors = []

for i, line6 in enumerate(lines6):
    contents = line6[:-1].split("\t\t")
    unixtime = float(contents[0])
    unixtime_err = float(contents[1])
    value = float(contents[2])
    value_err = float(contents[3])
    Xe131mUnixtimes.append(unixtime)
    Xe131mUnixtimeErrors.append(unixtime_err)
    Xe131mELifeValues.append(value)
    Xe131mELifeValueErrors.append(value_err)

#######################################
## Get the prediction lists
#######################################
# A simple method for correction the prediction
def SimpleCorrection(Elife, FitterUsed, Actual):
    TrueElife = pow( 1./Elife + 1./FitterUsed - 1./Actual, -1. )
    return TrueElife


fin2 = open(PredictionFile)
lines2 = fin2.readlines()
fin2.close()

UnixTimes2 = []
PredictedELifes = []
PredictedELifeLowers = []
PredictedELifeUppers = []
PredictedELifeLowerErrors = []
PredictedELifeUpperErrors = []
for i, line in enumerate(lines2):
    contents = line[:-1].split("\t\t")
    unixtime = float(contents[0])
    elife = float(contents[1])
    elife_lower = float(contents[2])
    elife_upper = float(contents[3])
    if unixtime>LastPointUnixTime+DaysAfterLastPoint*3600.*24.:
        break
    UnixTimes2.append(unixtime)
    PredictedELifes.append(elife)
    PredictedELifeLowers.append(elife_lower)
    PredictedELifeUppers.append(elife_upper)
    PredictedELifeLowerErrors.append( (elife_lower - elife)/elife)
    PredictedELifeUpperErrors.append((elife_upper - elife)/elife)
#    PredictedELifes.append( SimpleCorrection(elife, FitterUsedS1ExponentialConstant, S1ExponentialConstant) )
#    PredictedELifeLowers.append(SimpleCorrection(elife_lower, FitterUsedS1ExponentialConstant, S1ExponentialConstant) )
#    PredictedELifeUppers.append(SimpleCorrection(elife_upper, FitterUsedS1ExponentialConstant, S1ExponentialConstant) )
#    PredictedELifeLowerErrors.append( (elife_lower - elife)/elife)
#    PredictedELifeUpperErrors.append((elife_upper - elife)/elife)

from scipy.interpolate import interp1d

PredictionInterpolator = interp1d(UnixTimes2, PredictedELifes)

###################################
## Get the residual of the data points
###################################
ELifeValueDeviations = []
ELifeValueDeviationErrors = []
RnELifeValueDeviations = []
RnELifeValueDeviationErrors = []
for unixtime, elife, elife_err in zip(UnixTimes, ELifeValues, ELifeValueErrors):
    prediction = PredictionInterpolator(unixtime)
    residual = elife - prediction
    ELifeValueDeviations.append(residual / prediction *100.)
    ELifeValueDeviationErrors.append( elife_err / prediction * 100.)
for unixtime, elife,elife_err in zip(RnUnixtimes, RnELifeValues, RnELifeValueErrors):
    if unixtime > 1484761892:
        continue
    prediction = PredictionInterpolator(unixtime)
    residual = elife - prediction
    RnELifeValueDeviations.append(residual / prediction * 100.)
    RnELifeValueDeviationErrors.append( elife_err / prediction * 100.)

###################################
## Calculate the uncertainties in the first science run
###################################
ScienceRunStartUnixtime = 1479772800
ScienceRunEndUnixtime = 1484731512


###################################
## convert unixtimes to dates
###################################


Dates = [dt.datetime.fromtimestamp(ts) for ts in UnixTimes]
RnDates = [dt.datetime.fromtimestamp(ts) for ts in RnUnixtimes]
KrDates = [dt.datetime.fromtimestamp(ts) for ts in KrUnixtimes]
Xe129mDates = [dt.datetime.fromtimestamp(ts) for ts in Xe129mUnixtimes]
Xe131mDates = [dt.datetime.fromtimestamp(ts) for ts in Xe131mUnixtimes]
DateErrorLowers = []
DateErrorUppers = []
for ts, ts_err in zip(UnixTimes, UnixTimeErrors):
    date = dt.datetime.fromtimestamp(ts)
    date_err_lower = date - dt.datetime.fromtimestamp(ts - ts_err)
    date_err_upper = dt.datetime.fromtimestamp(ts + ts_err) - date
    DateErrorLowers.append( date_err_lower )
    DateErrorUppers.append( date_err_upper )
RnDateErrorLowers = []
RnDateErrorUppers = []
for ts, ts_err in zip(RnUnixtimes, RnUnixtimeErrors):
    date = dt.datetime.fromtimestamp(ts)
    date_err_lower = date - dt.datetime.fromtimestamp(ts - ts_err)
    date_err_upper = dt.datetime.fromtimestamp(ts + ts_err) - date
    RnDateErrorLowers.append( date_err_lower )
    RnDateErrorUppers.append( date_err_upper )
Xe129mDateErrorLowers = []
Xe129mDateErrorUppers = []
for ts, ts_err in zip(Xe129mUnixtimes, Xe129mUnixtimeErrors):
    date = dt.datetime.fromtimestamp(ts)
    date_err_lower = date - dt.datetime.fromtimestamp(ts - ts_err)
    date_err_upper = dt.datetime.fromtimestamp(ts + ts_err) - date
    Xe129mDateErrorLowers.append( date_err_lower )
    Xe129mDateErrorUppers.append( date_err_upper )
Xe131mDateErrorLowers = []
Xe131mDateErrorUppers = []
for ts, ts_err in zip(Xe131mUnixtimes, Xe131mUnixtimeErrors):
    date = dt.datetime.fromtimestamp(ts)
    date_err_lower = date - dt.datetime.fromtimestamp(ts - ts_err)
    date_err_upper = dt.datetime.fromtimestamp(ts + ts_err) - date
    Xe131mDateErrorLowers.append( date_err_lower )
    Xe131mDateErrorUppers.append( date_err_upper )
#UnixTimes2 = np.asarray(UnixTimes2)
#PredictedELifes = np.asarray(PredictedELifes)
#PredictedELifeLowers = np.asarray(PredictedELifeLowers)
#PredictedELifeUppers = np.asarray(PredictedELifeUppers)
#Dates2 = [dt.datetime.fromtimestamp(ts) for ts in UnixTimes2[UnixTimes2 < 1484731512]]
Dates2 = [dt.datetime.fromtimestamp(ts) for ts in UnixTimes2]
#UnixtimeOther = 1484731512 + 2.5*24*3600.
#Dates3 = [dt.datetime.fromtimestamp(ts) for ts in UnixTimes2[UnixTimes2 > UnixtimeOther]]


##############################
## Draw plot
##############################
from matplotlib import gridspec

XLimLow = dt.datetime.fromtimestamp(FirstPointUnixTime)
XLimLow = dt.datetime.fromtimestamp(1485802500)
XLimUp = dt.datetime.fromtimestamp(LastPointUnixTime+DaysAfterLastPoint*3600.*24.)


fig = plt.figure(figsize=(25.0, 16.0))
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

gs1 = gridspec.GridSpec(3,1)
ax = plt.subplot(gs1[0:3,:])

xfmt = md.DateFormatter('%Y-%m-%d')
ax.xaxis.set_major_formatter(xfmt)
ax.errorbar(Dates, ELifeValues, xerr=[DateErrorLowers,DateErrorUppers], yerr=[ELifeValueErrors,ELifeValueErrors], fmt='o', color='k', label="electron lifetime data points (S2/S1 method)")
ax.errorbar(RnDates, RnELifeValues,  xerr = [RnDateErrorLowers,RnDateErrorUppers], yerr=[RnELifeValueErrors,RnELifeValueErrors], fmt='o', color='deeppink', label="electron lifetime data points (from Rn analysis)")
ax.errorbar(KrDates, KrELifeValues, yerr = [KrELifeValueErrors, KrELifeValueErrors], fmt = 'o', color = 'g', label = "electron lifetime data points(from Kr83m analysis)")
ax.errorbar(Xe129mDates, Xe129mELifeValues,  xerr = [Xe129mDateErrorLowers,Xe129mDateErrorUppers], yerr=[Xe129mELifeValueErrors,Xe129mELifeValueErrors], fmt='o', color='darkmagenta', label="electron lifetime data points (from Xe129m analysis)")
ax.errorbar(Xe131mDates, Xe131mELifeValues,  xerr = [Xe131mDateErrorLowers,Xe131mDateErrorUppers], yerr=[Xe131mELifeValueErrors,Xe131mELifeValueErrors], fmt='o', color='darkorange', label="electron lifetime data points (from Xe131m analysis)")
ax.plot(
            Dates2,
            PredictedELifes,
            linewidth=2.,
            color = 'r',
            label='Best-fit trend',
           )
ax.fill_between(
                         Dates2,
                         PredictedELifeLowers,
                         PredictedELifeUppers,
                         color='b',
                         label='$\pm 1 \sigma$ C.L. region',
                         alpha=0.5,
                        )

#ax.plot(
#            Dates3,
#            PredictedELifes[UnixTimes2 > UnixtimeOther],
#            linewidth=2.,
#            color = 'r',
#            label='Best-fit trend',
#           )
#ax.fill_between(
#                         Dates3,
#                         PredictedELifeLowers[UnixTimes2 > UnixtimeOther],
#                         PredictedELifeUppers[UnixTimes2 > UnixtimeOther],
#                         color='b',
#                         label='$\pm 1 \sigma$ C.L. region',
#                         alpha=0.5,
#                        )


# plot the vertical lines for system change
ax.axvline( x=dt.datetime.fromtimestamp(1465937520), # first power outage
                    ymin = 0,
                    ymax = 650, 
                    linestyle = "--",
                    linewidth=3,
                    color='k',
                   )
ax.axvline( x=dt.datetime.fromtimestamp(1468597800), # LN2 test. PTR1 warm-up
                    ymin = 0,
                    ymax = 650, 
                    linestyle = "--",
                    linewidth=3,
                    color='k',
                   )
ax.axvline( x=dt.datetime.fromtimestamp(1484731512), # earthquake
                    ymin = 0,
                    ymax = 650, 
                    linestyle = "--",
                    linewidth=3,
                    color='k',
                   )



# fill the region
Xs = [
          dt.datetime.fromtimestamp(1471880000),
          dt.datetime.fromtimestamp(1472800000)
         ]
YLs = [0, 0]
YUs = [650, 650]
ax.fill_between(Xs, YLs, YUs, color='coral', alpha=0.7)
Xs = [
          dt.datetime.fromtimestamp(1475180000),
          dt.datetime.fromtimestamp(1475680000)
         ]
ax.fill_between(Xs, YLs, YUs, color='m', alpha=0.3)


Xs = [
          dt.datetime.fromtimestamp(ScienceRunStartUnixtime),
          dt.datetime.fromtimestamp(ScienceRunEndUnixtime)
         ]
YLs = [0, 0]
YUs = [650, 650]
ax.fill_between(Xs, YLs, YUs, color='coral', alpha=0.5)



# plot the text
ax.text( # Power outage
            dt.datetime.fromtimestamp(1465937520+2.*3600.*24.), 
            200., 
            'PTR warm-up',
            color='k',
            size=22.,
            rotation='vertical',
            )
ax.text( # LN2 test, PTR warm up
            dt.datetime.fromtimestamp(1468597800+2.*3600.*24.), 
            450., 
            'LN2 cooling test; PTR warm-up',
            color='k',
            size=22.,
            rotation='vertical',
            )
#ax.text( # Earthquake @ 01/18/2017
#            dt.datetime.fromtimestamp(1484731512+2.*3600.*24.), 
#            450., 
#            'Earthquake @ 01/18/17',
#            color='k',
#            size=22.,
#            rotation='vertical',
#            )
ax.text( # Gas-only circulation
            dt.datetime.fromtimestamp(1471880000-7.*3600.*24.), 
            675., 
            'Gas-only circulation',
            color='coral',
            size=22.,
            #rotation='vertical',
            )
ax.text(dt.datetime.fromtimestamp(1471880000), 580+20, "20 SLPM", color='coral', size=22.)
ax.text( # PUR upgrade
            dt.datetime.fromtimestamp(1475180000-5.*3600.*24.), 
            660., 
            'PUR upgrade',
            color='m',
            size=22.,
            alpha=0.5,
            #rotation='vertical',
            )


# text the flow rate
ax.text( dt.datetime.fromtimestamp(1464000000), 580+40., "$\sim$ 40 SLPM", size=20.,color='k')
ax.text( dt.datetime.fromtimestamp(1466500000), 580+20, "$\sim$ 55 SLPM", size=20.,color='k')
ax.text( dt.datetime.fromtimestamp(1469500000), 580+40, "45 - 50 SLPM", size=20.,color='k')
ax.text( dt.datetime.fromtimestamp(1473500000), 580+40, "$\sim$ 40 SLPM", size=20.,color='k')
ax.text( dt.datetime.fromtimestamp(1475700000), 580+20, "$\sim$ 54 SLPM", size=20.,color='k')

#ax.grid(True)
#XLimLow = datetime.datetime(2017, 3, 8, 0, 0)
#XLimUp = datetime.datetime(2017, 4, 28, 0, 0)
#XLimLow = datetime.datetime(2016, 11, 17, 0, 0)
#XLimUp = datetime.datetime(2017, 1, 22, 0, 0)
ax.set_xlim([XLimLow, XLimUp])
#ax.set_ylim([0, 650])
ax.set_ylim([400, 650])
#ax.set_ylim([300, 550])
ax.legend(loc = 'lower right',prop={'size':20})
ax.set_xlabel('Date', fontsize=30)
ax.set_ylabel('Electron lifetime $[\\mu s]$', fontsize=30)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)




ax.set_xlim([XLimLow, XLimUp])

fig.autofmt_xdate()

plt.savefig(FigureSaveName+".png", format='png')
plt.savefig(FigureSaveName+".pdf", format='pdf')

plt.show()
