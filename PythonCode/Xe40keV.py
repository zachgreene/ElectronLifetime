import numpy as np
import pandas as pd

import root_pandas
from root_pandas import read_root

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import matplotlib.dates as md

import hax

import lax
#from lax.lichens import sciencerun0
from lax.lichens import sciencerun1

from Tools import *
from FitElectronLifetime import FitElectronLifetime

hax.init(experiment='XENON1T')
dsets = hax.runs.datasets

XeIsotope = 'Xe40keV'

g1 = 0.1469
g2_b = 11.13
w = 0.0137
rectSelection = (150, 550, 4e3, 1.5e4)

#LifetimeFile =  '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWithXe40keV.txt'
LifetimeFile =  '/home/zgreene/xenon1t/ElectronLifetime/FitData/ElectronLifetimeDataWith' + XeIsotope + '.txt'
PathToData = '/project/lgrandi/zgreene/xenon1t/rootfiles/pax_v6.6.5/'
#PathToFigures = '/home/zgreene/xenon1t/ElectronLifetime/Figures/Xe40keV/'
PathToFigures = '/home/zgreene/xenon1t/ElectronLifetime/Figures/' + XeIsotope + '/'

bOverwrite = False


dSourceDates = dict(
                    Kr83m = [
                            [1489144526, 1489244895],
                            [1491207468, 1491300817],
#                            [1492418717, 1492514479]
                            ],
                    Cs137 = [
                            [1491312264, 1491342302]
                            ],
                    Th228 = [
                            [1491379392, 1491406193]
                            ],
                    Rn220 = [
                            [1489402529, 1489671382]
                            ],
                    AmBe = [
                            [1489671382, 1490960046]
                            ]
                    )

dSourceColors = dict(
                    Kr83m = 'g',
                    Cs137 = 'orangered',
                    Th228 = 'purple',
                    Rn220 = 'dodgerblue',
                    AmBe = 'r'
                    )


#cuts = [sciencerun0.S2SingleScatter()]
cuts = [sciencerun1.S2SingleScatter(),
        sciencerun1.S1SingleScatter(),
        sciencerun1.S2AreaFractionTop(),
        sciencerun1.S1AreaFractionTop(),
        sciencerun1.S2Width()]

cutnames = ['Cut' + cut.__class__.__name__ for cut in cuts]


AmBeFiles = [
                'BG_03_12_13',
                'BG_03_17_18',
                'BG_03_18_19',
                'BG_03_19_20',
                'BG_03_20_21',
                'BG_03_21_22',
                'BG_03_22_23',
                'BG_03_23_24',
                'BG_03_24_25',
                'BG_03_25_26',
                'BG_03_26_27',
                'BG_03_27_28',
                'BG_03_28_29',
                'BG_03_29_30',
                'BG_03_30_31',
                'BG_04_01_02',
                'BG_04_02_03',
                'BG_04_06_07',
                'BG_04_07_08',
                'BG_04_08_09',
                'BG_04_09_10',
                ]

NGFiles = [
                '2017_05_24_24hrs',
                '2017_05_25_24hrs',
                '2017_05_26_24hrs',
                '2017_05_27_24hrs',
                'BG_05_28_29',
                'BG_05_29_30',
                '2017_05_30_24hrs',
                ]

#AmBeFiles = ['AmBe']
NGFiles = ['NG']

NR_Type = 'AmBe'
#NR_Type = 'NG'

if NR_Type == 'AmBe':
    NR_Files = AmBeFiles
    print('\nLoading AmBe files\n')
elif NR_Type == 'NG':
    NR_Files = NGFiles
    print('\nLoading NG files\n')
else:
    NR_Files = AmBeFiles + NGFiles
    print('\nLoading AmBe + NG files\n')

PathToFigures = PathToFigures + NR_Type + '_SR1/'
if not os.path.exists(PathToFigures):
    os.mkdir(PathToFigures)

data = []
listOfFiles = []

################################
## load data
################################
#for i,File in enumerate(sorted(AmBeFiles)):
for i,File in enumerate(sorted(NR_Files)):
    data_temp = read_root(PathToData + File + '.root')
    NumEventsBeforeCuts = len(data_temp)
    data_temp = data_temp[(data_temp.s1 < 3800) & (data_temp.cs2_bottom < 2.2e5) &
                            (data_temp.r < 36.94)]

    print('Loaded %s with %i events - after cuts %i (%.2f%%) left'
            %(File + '.root', NumEventsBeforeCuts, len(data_temp),
            100.*len(data_temp)/NumEventsBeforeCuts)
            )
    listOfFiles.append(data_temp.run_number.unique())

    data.append(data_temp)

data = pd.concat(data)
data['dt'] = data['drift_time'] / 1000.
data['s2_bottom'] = data['s2'] * (1 - data['s2_area_fraction_top'])
data['E'] = (data['cs1'] / g1 + data['cs2_bottom'] / g2_b) * w

print(a)


#### merge with hax runs database to get start/end times ####
data = data.merge(dsets, left_on='run_number', right_on='number', how='left')


for cut in cuts:
    data = cut.process(data)

data_backup = data
#print(a)

#X_AmBe_FV = (97. - data.x)**2 + (43.5 - data.y)**2 + (-50. - data.z)**2 < 103.5**2
#X_AmBe_FV = (97. - data.x)**2 + (43.5 - data.y)**2 < 103.5**2

#X_NG_FV = (31.6 - data.x)**2 + (86.8 - data.y)**2 + (-50 - data.z)**2 < 111.5**2
#X_NG_FV = (31.6 - data.x)**2 + (86.8 - data.y)**2  < 111.5**2

#data = data[X_AmBe_FV]
#data = data[X_NG_FV]

#X1 = (data.r**2 < 36.94**2) & (data.CutS2AreaFractionTop) & (data.CutS2Width) & (data.CutS1AreaFractionTop) & (data.CutS2SingleScatter) & (data.CutS1SingleScatter)
#X1 = (data.r**2 < 36.94**2) & (data.CutS2AreaFractionTop) & (data.CutS2Width) & (data.CutS1AreaFractionTop) & (data.CutS1SingleScatter)
#X1 = (data.r**2 < 36.94**2) & (data.CutS2AreaFractionTop) & (data.CutS2Width) & (data.CutS1AreaFractionTop) & (data.CutS2SingleScatter)
X1 = (data.r**2 < 36.94**2) & (data.CutS2AreaFractionTop) & (data.CutS2Width) & (data.CutS1AreaFractionTop)
Xz =  (data.z > -96) & (data.z < -3)

X_NG = (data.source__type.str.lower() == 'neutron_generator')
X_AmBe = (data.source__type.str.lower() == 'ambe')
X_NR = (X_NG | X_AmBe) & X1


bFitEllipse = True
#fig, ax = plt.figure(figsize=(12, 6))


################################
## cs1 vs. cs2_bottom
################################
fig2, ax2 = plt.subplots(1, 1)
#ax2.hist2d(data.cs1[X_AmBe], data.cs2_bottom[X_AmBe], bins=(100, 100), range=((0, 2500), (0, 1e5)),
#            norm=colors.LogNorm(), cmap='viridis')
ax2.hist2d(data.cs1[X_NR], data.cs2_bottom[X_NR], bins=(100, 100), range=((0, 2500), (0, 1e5)),
            norm=colors.LogNorm(), cmap='viridis')
ax2.text(100, 1.7e4, r'$\rm{39.6 \,\, keV}$', rotation='horizontal',
            backgroundcolor='white', fontsize=12)
ax2.text(500, 2.5e4, r'$\rm{80.2 \,\, keV}$', rotation='horizontal',
            backgroundcolor='white', fontsize=12)
ax2.text(1300, 2.1e4, r'$\rm{163.9 \,\, keV}$', rotation='horizontal',
            backgroundcolor='white', fontsize=12)
ax2.text(1700, 3.5e4, r'$\rm{236.1 \,\, keV}$', rotation='horizontal',
            backgroundcolor='white', fontsize=12)

ax2.yaxis.get_major_formatter().set_powerlimits((1,2))
ax2.set_xlabel('$\\rm{cs1 \,\, [pe]}$', fontsize=20)
ax2.set_ylabel('$\\rm{cs2 \,\, bottom \,\, [pe]}$', fontsize=20)
#plt.savefig(PathToFigures + 'cs1_cs2_bottom_initial_' + NR_Type + '.png', format='png')
#plt.savefig(PathToFigures + 'cs1_cs2_bottom_initial_' + NR_Type + '.pdf', format='pdf')
#plt.show()


################################
## cs1 vs. cs2_bottom (Xe 40keV)
################################
fig3, ax3 = plt.subplots(1, 1)
ax3.hist2d(data.cs1[X_NR & Xz], data.cs2_bottom[X_NR & Xz], bins=(100, 100), range=((0, 1000), (0, 2.5e4)),
            norm=colors.LogNorm(), cmap='viridis')

if bFitEllipse:
    Bkg, Amp, MeanX, MeanY, Sigma1, Sigma2, Theta = FitEllipse(data[X_NR & Xz].cs1, data[X_NR & Xz].cs2_bottom, XeIsotope=XeIsotope)
    selection = Ellipse((MeanX, MeanY), 4.5*Sigma1, 4.5*Sigma2, Theta*180./np.pi, color='r', linewidth=3, fill=False)
    X = InsideEllipse(data.cs1, data.cs2_bottom, MeanX, MeanY, Sigma1, Sigma2, Theta) < 4.5
    NR_Selection = 'ellipse'
else:
    cs1Min = rectSelection[0]
    cs1Max = rectSelection[1]
    cs1Width = cs1Max - cs1Min
    cs2Min = rectSelection[2]
    cs2Max = rectSelection[3]
    cs2Width = cs2Max - cs2Min
    selection = Rectangle((cs1Min, cs2Min), cs1Width, cs2Width, linewidth=3,
                            fill=False, facecolor='none', edgecolor='orangered')
    X = (data.cs1 > cs1Min) & (data.cs1 < cs1Max) & (data.cs2_bottom > cs2Min) & (data.cs2_bottom < cs2Max)
    NR_Selection = 'rectangle'

ax3.add_artist(selection)
ax3.set_xlabel('$\\rm{cs1 \,\, [pe]}$', fontsize=20)
ax3.set_ylabel('$\\rm{cs2 \,\, bottom \,\, [pe]}$', fontsize=20)
#plt.savefig(PathToFigures + 'cs1_cs2_bottom_initial_40kev_' + NR_Type + '_' + NR_Selection + '.png', format='png')
#plt.savefig(PathToFigures + 'cs1_cs2_bottom_initial_40kev_' + NR_Type + '_' + NR_Selection + '.pdf', format='pdf')
plt.show()

X2 = X & X1


fig, ax = plt.subplots(1, 1, figsize=(10, 6))
plt.hist2d(data.dt[X_NR], data.cs1[X_NR], bins=(120, 120), range=((0, 730), (0, 2500)),
            norm=colors.LogNorm(), cmap='viridis', alpha=0.5)
cb = plt.colorbar()
cb.remove()
vmin, vmax = plt.gci().get_clim()
plt.hist2d(data.dt[X], data.cs1[X], bins=(120, 120), range=((0, 730), (0, 2500)),
            norm=colors.LogNorm(), cmap='viridis', vmin=vmin, vmax=vmax)
ax.text(735, 310, '39.6 keV', rotation='horizontal')
ax.text(735, 670, '80.2 keV', rotation='horizontal')
ax.text(735, 1190, '163.9 keV', rotation='horizontal')
ax.text(735, 1700, '236.1 keV', rotation='horizontal')
ax.set_xlabel('$\\rm{drift \,\, time \,\, [\mu s]}$', fontsize=20)
ax.set_ylabel('$\\rm{cs1 \,\, [pe]}$', fontsize=20)
plt.subplots_adjust(right=0.9)
#plt.savefig(PathToFigures + 'dt_cs1_' + NR_Type + '_' + NR_Selection + '.png', format='png')
#plt.savefig(PathToFigures + 'dt_cs1_' + NR_Type + '_' + NR_Selection + '.pdf', format='pdf')
#plt.show()


#data = data_keep
data_ces = data[X_NR & (data.dt > 60) & (data.dt < 710)]
Energies = [39.6, 80.2, 163.9, 236.1]

Isotopes = ['$\\rm{^{129}Xe\,\,\,\,39.4\,\,keV}$',
            '$\\rm{^{131}Xe\,\,\,\,80.2\,\,keV}$',
            '$\\rm{^{131m}Xe\,\,\,\,163.9\,\,keV}$',
            '$\\rm{^{129m}Xe\,\,\,\,236.1\,\,keV}$']

fig5, ax5 = plt.subplots(1, 1, figsize=(10, 6))
hist, bin_edges = np.histogram(data_ces.E[data_ces.cs2_bottom > 4e3], bins=np.linspace(0, 250, 251))
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
bin_width = bin_centers[1] - bin_centers[0]
plt.errorbar(bin_centers, hist, xerr=bin_width, yerr=np.sqrt(hist),
                color='darkgreen', fmt='o', markeredgewidth=0, markersize=2)

for Isotope,Energy in zip(Isotopes,Energies):
    plt.axvline(x=Energy, color='r')
#    plt.text(Energy+2, 6e3, '%s %.1f keV' %(Isotope, Energy), rotation='vertical', color='r')
#    print(Energy)
    plt.text(Energy-8, 7e3, '%s' %Isotope, rotation='vertical', color='r', fontsize=14)

plt.gca().set_xlim(0, 250)
plt.gca().set_yscale('log')
plt.gca().set_xlabel('energy [keV]')
plt.gca().set_ylabel('counts')
plt.tight_layout()
#plt.savefig(PathToFigures + 'ambe_ces_' + NR_Type + '_' + NR_Selection + '.png', format='png')
#plt.savefig(PathToFigures + 'ambe_ces_' + NR_Type + '_' + NR_Selection + '.pdf', format='pdf')
#plt.show()

#AmBeLists = []
#NumAmBeLists = 2
NRLists = []
NumNRLists = 2
data_all_time = data[X2 & Xz]
# get events/day
livetimes = []
livetime_errs = []
counts = []
count_errs =[]
datetimes = []
datetime_low_errs = []
datetime_up_errs = []

Dates = pd.to_datetime(
                        np.unique(
                                    pd.DatetimeIndex(data.start.unique()).strftime('%Y-%m-%d')
                                    )
                        )

for Date in Dates:
    DatetimeStart = Date
    DatetimeEnd = Date + datetime.timedelta(days=1)
    subset = dsets[(dsets.start >= DatetimeStart) & (dsets.end < DatetimeEnd) &
                    (dsets.number.isin(data_backup.run_number))]
    del_time = subset.end - subset.start
    livetime = del_time.sum().total_seconds()
    count = len(data_all_time[(data_all_time.start >= DatetimeStart) & (data_all_time.start < DatetimeEnd)])
    livetimes.append(livetime)
    counts.append(count / (livetime / 3600. / 24.))
    count_errs.append(count**0.5 / (livetime / 3600. / 24.))
    dtime = subset.start.min() + del_time.sum() / 2.
    datetimes.append(dtime)
    datetime_low_errs.append(del_time.sum() / 2)
    datetime_up_errs.append(subset.end.max() - dtime)
    subset = subset[(subset.source__type.str.lower() == 'ambe') |
                    (subset.source__type.str.lower() == 'neutron_generator')]
    if (('ambe' in subset.source__type.str.lower().unique()) or
                            ('neutron_generator' in subset.source__type.str.lower().unique())):
        NRLists.append(subset.number.values)


#for i,Fileset in enumerate(listOfFiles):
#    subset = dsets[dsets.number.isin(Fileset)]
#    print(subset.number.values)
#    del_time = subset.end - subset.start
#    livetime = del_time.sum().total_seconds()
#    count = len(data_all_time[data_all_time.run_number.isin(Fileset)])
#    livetimes.append(livetime)
#    counts.append(count / (livetime / 3600. / 24.))
#    count_errs.append(count**0.5 / (livetime / 3600. / 24.))
#    dtime = subset.start.min() + del_time.sum() / 2.
#    datetimes.append(dtime)
#    datetime_low_errs.append(del_time.sum() / 2)
#    datetime_up_errs.append(subset.end.max() - dtime)
#    if 'ambe' in subset.source__type.str.lower().unique():
#        AmBeLists.append(Fileset)

Colors = ['darkmagenta', 'darkorange']

#plt.gca().set_yscale('log', nonposy='clip')
fig = plt.figure(figsize=(18,8))
xfmt = md.DateFormatter('%Y-%m-%d')
plt.gca().xaxis.set_major_formatter(xfmt)

plt.errorbar(datetimes, counts, xerr=[datetime_low_errs, datetime_up_errs], yerr=count_errs,
#                color='darkgreen', fmt='o', markeredgewidth=0, label='$\\rm{^{129}Xe\,\, 39.6\,\,keV}$')
                color='darkgreen', fmt='o', markeredgewidth=0, label='selected events')
yLims = plt.ylim()
TextHeight = np.mean(yLims)
if True:
    for key,value in dSourceDates.items():
        for unixtimes in value:
                plt.gca().fill_between([datetime.datetime.fromtimestamp(unixtimes[0]),
                                        datetime.datetime.fromtimestamp(unixtimes[1])],
                                        [0, 0], [yLims[1], yLims[1]],
                                        color=dSourceColors[key], alpha=0.3)
                plt.gca().text(datetime.datetime.fromtimestamp(unixtimes[0] + 3*3600),
                                TextHeight, key, color=dSourceColors[key],
                                rotation='vertical')

plt.xlabel('Date')
plt.ylabel('Rate [day$^{-1}$]')
plt.legend(loc='upper right', prop={'size': 12})
days = md.DayLocator()
plt.gca().xaxis.set_minor_locator(days)
fig.autofmt_xdate()
plt.tight_layout()
#plt.savefig(PathToFigures + 'daily_rate_40kev_' + NR_Type + '_' + NR_Selection + '.png', format='png')
#plt.savefig(PathToFigures + 'daily_rate_40kev_' + NR_Type + '_' + NR_Selection + '.pdf', format='pdf')
#plt.show()



#LengthOfFilelists = int(len(AmBeLists) / NumAmBeLists)
#AmBeListsToFit = []
LengthOfFilelists = int(len(NRLists) / NumNRLists)
NRListsToFit = []
if 'NG' in NR_Type:
    NRListsToFit.append(np.concatenate(NRLists[:4]))
    NRListsToFit.append(np.concatenate(NRLists[4:]))
elif 'AmBe' in NR_Type:
    for i in range(NumNRLists):
        IndexBegin = i*LengthOfFilelists
        IndexEnd = (i+1)*LengthOfFilelists
        if IndexEnd >= len(NRLists):
            NRListsToFit.append(np.concatenate(NRLists[IndexBegin:None]))
        else:
            NRListsToFit.append(np.concatenate(NRLists[IndexBegin:IndexEnd]))


data_keep = data
X3 = X & X_NR
#data = data[X3]


#plt.show()
#print(a)
plt.close(plt.gcf())

print(a)

NumEventsTot = 0
NumEventsFit = 0

#for i in range(len(NRListsToFit):
fig_lifetime, ax_lifetime = plt.subplots(figsize=(14, 8))
data = data_keep[X3].loc[data_keep[X3].number.isin(NRListsToFit[i])]
(
    ax_lifetime,
    MeanUnixtime,
    ErrUnixtime,
    Lifetime,
    LifetimeErr,
    NumEvents) = FitElectronLifetime(data,
                                        's2_bottom',
                                        ax_lifetime,
                                        FitEstimates=[1e4, 500]
                                        )
NumEventsTot += len(data)
NumEventsFit += NumEvents
print(NumEventsTot, NumEventsFit)

MeanDatetime = datetime.datetime.fromtimestamp(MeanUnixtime)

#plt.savefig(PathToFigures + 'lifetime_' + NR_Type + '_' + NR_Selection + '_' +
#            MeanDatetime.strftime('%y%m%d') + '.png', format='png')
#plt.savefig(PathToFigures + 'lifetime_' + NR_Type + '_' + NR_Selection + '_' +
#            MeanDatetime.strftime('%y%m%d') + '.pdf', format='pdf')
plt.show()

print(a)

#if bOverwrite:
#ValuesToInsert = (MeanUnixtime, ErrUnixtime, Lifetime, LifetimeErr)
#InsertLifetimeToFile(ValuesToInsert, LifetimeFile,  bOverwrite=bOverwrite)

for i in range(len(NRListsToFit)):
#    data = data[X3]
    data = data_keep[X3].loc[data_keep[X3].number.isin(NRListsToFit[i])]
    x = data.dt
    y = data.s2_bottom
    #x = data[X].dt
    #y = data[X].s2_bottom
    #data = data[X]
    DriftTimeTPC = 727.
    MinDriftTime = 40
    MaxDriftTime = 710
    NumBinsDriftTime = 30
    BinsOmitLow = 0
    BinsOmitUp = 0
    FactorAbove = 1.4
    FactorBelow = 0.6
    qs = [16, 50, 84]
    #xBins = np.linspace(0.8*MinDriftTime, 1.1*MaxDriftTime, NumBinsDriftTime)
    xBins = np.linspace(MinDriftTime, MaxDriftTime, NumBinsDriftTime)
    #yBins = ParameterManager.GetBinsForPrelimFit(S2ToFit)
    yBins = np.linspace(1000, 15000, 100)
    # get preliminary fit
    try:
        PrelimOpts = GetPreliminaryFit(x, y, xBins, qs=qs, FitType='Exp', BinsOmitLow=BinsOmitLow, BinsOmitUp=BinsOmitUp)
    except:
        print('error with prelim 100 bins, trying with 80')
        PrelimOpts = GetPreliminaryFit(x, y, np.linspace(0, MaxDriftTime, 80), qs=qs, FitType='Exp', BinsOmitLow=BinsOmitLow, BinsOmitUp=BinsOmitUp)
        print('prelim fit worked with 80 bins')
    p0, p1 = PrelimOpts
    print(p0, p1)
    # exclude data far outside of preliminary fit
    dataOmit = data[
                    (data.s2_bottom < FactorBelow * p0 * np.exp(-data.dt / p1)) |
                    (data.s2_bottom > FactorAbove * p0 * np.exp(-data.dt / p1))
                    ]
    data = data[data.s2_bottom > FactorBelow * p0 * np.exp(-data.dt / p1)]
    data = data[data.s2_bottom < FactorAbove * p0 * np.exp(-data.dt / p1)]
    TimePoints = np.linspace(0, DriftTimeTPC, DriftTimeTPC+1)
    # to draw as s2-dt boundaries on plot
    PrelimFitLine = FitExponential(TimePoints, *PrelimOpts)
    (BinCenters,
    BinErrs,
    OptParams,
    OptCovs) = GetBinGaussianParameters(
                        data.dt,
                        data.s2_bottom,
                        xBins=np.linspace(MinDriftTime, MaxDriftTime, NumBinsDriftTime),
                        sPathToFileSave='/home/zgreene/fig.png',
                        bUseROOT=True
                        )
    Means = OptParams.T[1]
    AvgErrs = np.sqrt(OptCovs.T[1][1])
    UpErrs = Means + AvgErrs
    LowErrs = Means - AvgErrs
    if BinsOmitLow == 0:
        bStatusFail, popt, pcov = FitExponentialROOT(
                                                        BinCenters[BinsOmitUp:],
                                                        Means[BinsOmitUp:],
                                                        BinErrs[BinsOmitUp:],
                                                        AvgErrs[BinsOmitUp:],
                                                        [1e4, 600]
                                                        )
    else:
        bStatusFail, popt, pcov = FitExponentialROOT(
                                                        BinCenters[BinsOmitUp:-BinsOmitLow],
                                                        Means[BinsOmitUp:-BinsOmitLow],
                                                        BinErrs[BinsOmitUp:-BinsOmitLow],
                                                        AvgErrs[BinsOmitUp:-BinsOmitLow],
                                                        [1e4, 600]
                                                        )
    Lifetime = popt[1]
    LifetimeErr = (pcov[1][1])**0.5
    S2Min = 1000
    S2Max = 15000
    fig_lifetime, ax_lifetime = plt.subplots(figsize=(14, 8))
    ax_lifetime = GetLifetimeLogPlot(data, dataOmit, DriftTimeTPC, ax_lifetime, 's2_bottom', S2Min=S2Min, S2Max=S2Max)
    ############################################
    ########## draw electon lifetime  ##########
    ############################################
    # draw s2-dt cut lines
    ax_lifetime.plot(TimePoints, FactorAbove*PrelimFitLine,
                linestyle='dashed', color='purple', label='Remove events above distribution')
    ax_lifetime.plot(TimePoints, FactorBelow*PrelimFitLine,
                linestyle='dashed', color='blue', label='Remove events below distribution')
    if BinsOmitLow == 0:
        ax_lifetime.errorbar(
                                BinCenters[BinsOmitUp:],
                                Means[BinsOmitUp:],
                                xerr=BinErrs[BinsOmitUp:],
                                yerr=AvgErrs[BinsOmitUp:],
                                fmt='o', color='orangered', capsize=0,
                                markeredgewidth=0, linewidth=2
                                )
    else:
        ax_lifetime.errorbar(
                                BinCenters[BinsOmitUp:-BinsOmitLow],
                                Means[BinsOmitUp:-BinsOmitLow],
                                xerr=BinErrs[BinsOmitUp:-BinsOmitLow],
                                yerr=AvgErrs[BinsOmitUp:-BinsOmitLow],
                                fmt='o', color='orangered', capsize=0,
                                markeredgewidth=0, linewidth=2
                                )
        # draw best fit lifetime
    FitDriftTimes = np.linspace(MinDriftTime, MaxDriftTime, 1000)
    ax_lifetime.step(FitDriftTimes, FitExponential(FitDriftTimes, *popt), linewidth=2, color='orangered', linestyle='-')
    
    #ax_lifetime.set_xlabel(r'$\rm{dt [}\mu \rm{s]}$', fontsize=20)
    ax_lifetime.set_ylabel('$\\rm{s2\,\,bottom\,\, [pe]}$', fontsize=20)
    ax_lifetime.set_yscale('log', fontsize=20)
    ax_lifetime.set_xlim(0, DriftTimeTPC)
    ax_lifetime.set_ylim(S2Min, S2Max)
    LifetimeText = 'Electron Lifetime = %.2f $\pm$ %.2f $\mu s$' % (Lifetime, LifetimeErr)

    ErrDatetime = (data.end.max() - data.start.min()) / 2.
    MeanDatetime = data.start.min() + ErrDatetime
    MeanUnixtime = time.mktime(MeanDatetime.timetuple())
    ErrUnixtime = ErrDatetime.total_seconds()

#    ValuesToInsert = (MeanUnixtime, ErrUnixtime, Lifetime, LifetimeErr)
#    InsertLifetimeToFile(ValuesToInsert, LifetimeFile,  bOverwrite=bOverwrite)






    NumEventsTot += len(data[(data.dt > MinDriftTime) & (data.dt < MaxDriftTime)])





    FigTitle = '%s $\pm$ %.1f days' %(MeanDatetime.strftime('%Y-%m-%d %H:%M'), ErrDatetime.days + ErrDatetime.seconds/3600. / 24.)
    
    ax_lifetime.text(40, 1.4*ax_lifetime.get_ylim()[0], LifetimeText, color='orangered',
                    fontsize=20, weight='bold', alpha=1)
#    ax_lifetime.text(40, 1.25*ax_lifetime.get_ylim()[0], '@ %.3f $\pm$ %.3f' % (MeanUnixtime, ErrUnixtime),
#                    color='orangered', fontsize=20, weight='bold', alpha=1)
    
#    ax_lifetime.legend(loc='best', fontsize=20)
    ax_lifetime.set_xlabel(r'$\rm{drift\,\, time\,\, [}\mu \rm{s]}$', fontsize=20)
    ax_lifetime.set_title(FigTitle, fontsize=20)
    plt.tight_layout()
    print('%.1f +/- %.1f' %(Lifetime, LifetimeErr))
#    plt.savefig(PathToFigures + 'lifetime_' + NR_Type + '_' + MeanDatetime.strftime('%y%m%d') + '.png', format='png')
#    plt.savefig(PathToFigures + 'lifetime_' + NR_Type + '_' + MeanDatetime.strftime('%y%m%d') + '.pdf', format='pdf')
    plt.show()

print('Length of data: %i' %len(data_keep[X3]))
print('NumEvents: %i' %NumEventsTot)
