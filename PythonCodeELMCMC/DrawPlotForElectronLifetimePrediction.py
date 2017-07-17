import numpy as np
import scipy as sp
from scipy.interpolate import interp1d

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib import gridspec
import datetime as dt
import time
import pickle
import sys
from Tools import *
import FormPars

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

Xe129mELifeDataFile = '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithXe129m.txt'
Xe131mELifeDataFile = '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithXe131m.txt'


#######################################
### Get single scatter elife data
#######################################
UnixTimes, UnixTimeErrors, ELifeValues, ELifeValueErrors = LoadFitData('SingleScatter', PathToFile=ELifeDataFile)

FirstPointUnixTime = UnixTimes[0]
LastPointUnixtime = UnixTimes[len(UnixTimes)-1]

######################################
## Get Rn elife data
######################################
RnUnixtimes, RnUnixtimeErrors, RnELifeValues, RnELifeValueErrors = LoadFitData('Rn', PathToFile=RnELifeDataFile)

LastPointUnixtime = RnUnixtimes[-1]

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
KrUnixtimes, KrUnixtimeErrors, KrELifeValues, KrELifeValueErrors = LoadFitData('Kr83', PathToFile=Kr83ELifeDataFile)

####################################
# Get Xe129m/Xe131m elifes 
####################################
Xe129mUnixtimes, Xe129mUnixtimeErrors, Xe129mELifeValues, Xe129mELifeValueErrors = LoadFitData('Xe129', PathToFile=Xe129mELifeDataFile)
Xe131mUnixtimes, Xe131mUnixtimeErrors, Xe131mELifeValues, Xe131mELifeValueErrors = LoadFitData('Xe131', PathToFile=Xe131mELifeDataFile)

#######################################
## Get the prediction lists
#######################################
(PredictionUnixtimes,
PredictedELifes,
PredictedELifeLowers,
PredictedELifeUppers,
PredictedELifeLowerErrors,
PredictedELifeUpperErrors) = LoadPredictions(PredictionFile,
                                            LastPointUnixtime,
                                            DaysAfterLastPoint=DaysAfterLastPoint)


PredictionInterpolator = interp1d(PredictionUnixtimes, PredictedELifes)

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
Dates2 = [dt.datetime.fromtimestamp(ts) for ts in PredictionUnixtimes]
#UnixtimeOther = 1484731512 + 2.5*24*3600.


##############################
## Draw plot
##############################
XLimLow = dt.datetime.fromtimestamp(FirstPointUnixTime)
#XLimLow = dt.datetime.fromtimestamp(1485802500)
XLimUp = dt.datetime.fromtimestamp(LastPointUnixtime+DaysAfterLastPoint*3600.*24.)


fig = plt.figure(figsize=(25.0, 16.0))
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

gs1 = gridspec.GridSpec(3,1)
ax = plt.subplot(gs1[0:3,:])

xfmt = md.DateFormatter('%Y-%m-%d')
ax.xaxis.set_major_formatter(xfmt)

ax.errorbar(Dates, ELifeValues, xerr=[DateErrorLowers,DateErrorUppers],
            yerr=[ELifeValueErrors,ELifeValueErrors], fmt='o', color='k',
            label="electron lifetime data points (S2/S1 method)")

ax.errorbar(RnDates, RnELifeValues,  xerr = [RnDateErrorLowers,RnDateErrorUppers],
            yerr=[RnELifeValueErrors,RnELifeValueErrors], fmt='o', color='deeppink',
            label="electron lifetime data points (from Rn analysis)")

ax.errorbar(KrDates, KrELifeValues, yerr = [KrELifeValueErrors, KrELifeValueErrors],
            fmt = 'o', color = 'g', label = "electron lifetime data points(from Kr83m analysis)")

ax.errorbar(Xe129mDates, Xe129mELifeValues, xerr=[Xe129mDateErrorLowers,Xe129mDateErrorUppers],
            yerr=[Xe129mELifeValueErrors, Xe129mELifeValueErrors], fmt='o', color='darkmagenta',
            label="electron lifetime data points (from Xe129m analysis)")

ax.errorbar(Xe131mDates, Xe131mELifeValues, xerr=[Xe131mDateErrorLowers,Xe131mDateErrorUppers],
            yerr=[Xe131mELifeValueErrors, Xe131mELifeValueErrors], fmt='o', color='darkorange',
            label="electron lifetime data points (from Xe131m analysis)")


CathodeVoltages = FormPars.GetCathodeVoltages()

# plot times when voltage is not 0 kV, otherwise fill
for CathodeVoltage in CathodeVoltages:
    Dates2 = [dt.datetime.fromtimestamp(ts) for ts in PredictionUnixtimes
                if(ts >= CathodeVoltage[0][0] and ts < CathodeVoltage[0][1])]

    if CathodeVoltage[1][0] == 0 and CathodeVoltage[1][1] == 0:
        ax.fill_between(Dates2, 0, 650, color='y', alpha=0.5, label=r'$V_{C} = 0$ kV')
        continue

    ELifesToPlot = [ELife for ts,ELife in zip(PredictionUnixtimes,PredictedELifes)
                    if (ts >= CathodeVoltage[0][0] and ts < CathodeVoltage[0][1])]
    ELifesLowToPlot = [ELife for ts,ELife in zip(PredictionUnixtimes,PredictedELifeLowers)
                    if (ts >= CathodeVoltage[0][0] and ts < CathodeVoltage[0][1])]
    ELifesUpToPlot = [ELife for ts,ELife in zip(PredictionUnixtimes,PredictedELifeUppers)
                    if (ts >= CathodeVoltage[0][0] and ts < CathodeVoltage[0][1])]
    ax.plot(
                Dates2,
#                PredictedELifes,
                ELifesToPlot,
                linewidth=2.,
                color = 'r',
                label='Best-fit trend',
               )
    ax.fill_between(
                             Dates2,
                             ELifesLowToPlot,
                             ELifesUpToPlot,
#                             PredictedELifeLowers,
#                             PredictedELifeUppers,
                             color='b',
                             label=r'$\pm 1 \sigma$ C.L. region',
                             alpha=0.5,
                            )



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
ax.set_ylim([0, 650])
#ax.set_ylim([400, 650])
#ax.set_ylim([300, 550])
from collections import OrderedDict

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys(), loc = 'lower right',prop={'size':20})
#ax.legend(loc = 'lower right',prop={'size':20})
ax.set_xlabel('Date', fontsize=30)
ax.set_ylabel('Electron lifetime $[\\mu s]$', fontsize=30)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)




ax.set_xlim([XLimLow, XLimUp])

fig.autofmt_xdate()

plt.savefig(FigureSaveName+".png", format='png')
plt.savefig(FigureSaveName+".pdf", format='pdf')

plt.show()
