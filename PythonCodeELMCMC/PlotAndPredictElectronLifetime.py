import ElectronLifetimeTrend
from ElectronLifetimeTrend import *

#import MyHistorianLib
#from MyHistorianLib import GetUnixTimeFromTimeStamp

import Tools
from Tools import *

import FormPars

import ROOT
from ROOT import TCanvas 
from ROOT import TFile
from ROOT import TTree
from ROOT import TChain
from ROOT import TCut
from ROOT import TH1
from ROOT import TH1D
from ROOT import TH2
from ROOT import TH2D
from ROOT import TGraph
from ROOT import TGraphErrors
from ROOT import TF1
from ROOT import TLatex
from ROOT import TLine
from ROOT import TLegend
from ROOT import TPad

import numpy as np
from scipy.optimize import minimize
from numpy.linalg import inv

if len(sys.argv)<2:
    print("======= Synatex ========")
    print("python PlotAndPredictElectronLifetime.py <Historian file>  <Fit result txt> <E-Life data txt> <Figure output filename> <Prediction txt output> <burn-in iteraction cut>")
    print("NOTICE: do remember to add \"-r\" at the end")
    exit()

HistorianFile = sys.argv[1]
FitResultInput = sys.argv[2]
#ElectronLifetimeDataFile = sys.argv[3]
#FigureOutputFile = sys.argv[4]
PredictionOutputFile = sys.argv[3]
BurnInCutOff = int(sys.argv[4])


NumOfInterpretation = 1000
MaximumELifePlot = 500.
MinimumELifePlot = 1.
NumOfTrials = 1000


S1ExponentialConstant = 2040 # us. 
# setting the parameters
#MinUnixTime = GetUnixTimeFromTimeStamp("05/17/16 00:00:00 ")
#MaxUnixTime = GetUnixTimeFromTimeStamp("05/31/17 00:00:00 ")
#MinUnixTime = GetUnixTimeFromTimeStamp("05/17/16 00:00:00 ")
#MaxUnixTime = GetUnixTimeFromTimeStamp("02/21/17 00:00:00 ")
MinUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMinTimeStamp())
MaxUnixTime = GetUnixTimeFromTimeStamp(FormPars.GetMaxTimeStamp()) + 50*24*3600
default_pars = FormPars.GetDefaultPars()

def RegulatePars(InputPars):
    OutputPars = []
    for par in InputPars:
        if par<=0:
            OutputPars.append(0)
        else:
            OutputPars.append(par)
    return OutputPars





# Need to load the electron lifetime measurement

#f1 = open(ElectronLifetimeDataFile)
#lines1 = f1.readlines()
#f1.close()

#ElectronLifetimeData = {}
#gEL = TGraphErrors()
##############################
## fill in the data
##############################
## electron lifetime
#UnixTimes = []
#Values = []
#ValueErrors = []
#for i, line in enumerate(lines1):
#    contents = line[:-1].split("\t\t")
#    unixtime = float(contents[0])
#    unixtime_err = float(contents[1])
#    value = float(contents[2])
#    value_err = float(contents[3])
#    # correct the e-life by the S1 term
#    value = value*S1ExponentialConstant / (S1ExponentialConstant - value)
#    value_err = value*np.sqrt( np.power(value_err/value,2.0)+np.power(value_err/ (S1ExponentialConstant - value), 2.) )
#    UnixTimes.append(unixtime)
#    Values.append(value)
#    ValueErrors.append(value_err)
#    gEL.SetPoint(i, unixtime, value)
#    gEL.SetPointError(i, unixtime_err, value_err)
#ElectronLifetimeData['UnixTimes'] = UnixTimes
#ElectronLifetimeData['Values'] = Values
#ElectronLifetimeData['ValueErrors'] = ValueErrors
#LastPointUnixTime = UnixTimes[len(UnixTimes)-1]



#############################
## Get the MCMC result
#############################
import pickle
MCMCResults = pickle.load(open(FitResultInput, 'rb'))
sampler = MCMCResults['sampler']
#ndim = MCMCResults['ndim']
#nwalkers = MCMCResults['nwalkers']
#niterations = MCMCResults['niterations']
ndim = sampler.__dict__['dim']
niterations = sampler.__dict__['iterations']
nwalkers = sampler.__dict__['k']
print(ndim, nwalkers, niterations)

# Cutting Burn-In and unreasonable region
# the value shall already be from "PlotElectronLifetime.py"
#samples_cut = GetBurnInCutoffSamples(MCMCResults['chain'], int(BurnInCutOff))
#samples_cut = GetBurnInCutoffSamples(sampler.__dict__['_chain'], int(BurnInCutOff))
samples_cut = GetBurnInCutoffSamplesV2(sampler.__dict__['_chain'], int(BurnInCutOff), int(BurnInCutOff)+199)

#####################################
## Calculate the best fit from MCMC results
#####################################
mean = np.average(samples_cut, axis=0)
print(mean)


#############################
# plot 
#############################
Pars, IfSth=FormPars.FormPars(mean)
#MaxUnixTime = LastPointUnixTime + 60.*3600.*24. # 2 month after the last data point
# The main Light yield Trend
pElectronLifetimeTrend = MyElectronLifetimeTrend(HistorianFile, MinUnixTime, MaxUnixTime, Pars)

# get the graphs 
UnixTimes = np.linspace(MinUnixTime, MaxUnixTime, NumOfInterpretation)
# include times of discontinuity
ImpactfulUnixtimes = FormPars.GetImpactfulUnixtimes()
for ImpactfulUnixtime in ImpactfulUnixtimes:
	if ImpactfulUnixtime-1 not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime-1]), ImpactfulUnixtime-1)
	if ImpactfulUnixtime not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime]), ImpactfulUnixtime)
	if ImpactfulUnixtime+1 not in UnixTimes:
		UnixTimes = np.insert(UnixTimes, len(UnixTimes[UnixTimes < ImpactfulUnixtime+1]), ImpactfulUnixtime+1)

#Taus = []
#gI = TGraph()
#gTau = TGraph()
#gFlow = TGraph()
#for i, unixtime in enumerate(UnixTimes):
#    I = pElectronLifetimeTrend.GetImpurityTrend().GetConcentrations(unixtime)[1]
#    #print(I)
#    Tau = pElectronLifetimeTrend.GetElectronLifetime(unixtime)
#    Taus.append(Tau)
#    LiquidFlow, GasFlow, CoolingPower, _ = pElectronLifetimeTrend.GetImpurityTrend().GetHistorianData().GetHistorian(unixtime)
#    gI.SetPoint(i, unixtime, I)
#    gTau.SetPoint(i, unixtime, Tau)
#    gFlow.SetPoint(i, unixtime, LiquidFlow)
    





#hframe1 = GetFrame("hframe1", MinUnixTime*3600*24., MaxUnixTime*3600*24., MinimumELifePlot, MaximumELifePlot, "Date", "Electron lifetime[#mus]")
#hframe1.GetXaxis().SetTimeDisplay(1)
#hframe1.GetXaxis().SetTimeFormat("%m/%d %F 1970-01-01 00:00:00")


#tex = TLatex()
#tex.SetTextFont(132)
#tex.SetTextSize(0.1)






#############################
# calculate:
# 1) Maximum
# 2) Reaching date for 90% capacity
# 3) Reaching date for 500us
#############################
# The electron lifetime Trend for varied parameters

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

# Initial the 1 sigma lower/upper of the trend
Trends = []
for i in range(len(UnixTimes)):
    Trends.append([])


# pickup randomly the pars in the MCMC walker
Pars_Trials = PickupMCMCPosteriors(samples_cut, NumOfTrials)
for i, pars_random in enumerate(Pars_Trials):
    print("i = "+str(i))
    pars, IfSth = FormPars.FormPars(pars_random)
    pElectronLifetimeTrend.SetParameters(pars)
    for j, unixtime in enumerate(UnixTimes):
        Trends[j].append(pElectronLifetimeTrend.GetElectronLifetime(unixtime))
StandardDeviation1Sigma = [15.4, 50., 84.6]
Taus = []
LowerBoundaries = []
UpperBoundaries = []
for i, deviation in enumerate(StandardDeviation1Sigma):
    for Values in Trends:
        Boundary = PercentileWithInf(Values, deviation)
        if i==0:
            LowerBoundaries.append(Boundary)
        elif i==1:
            Taus.append(Boundary)
        else:
            UpperBoundaries.append(Boundary)
#print(MaximumELifes)




######### moved from bottom because terminal cannot show root plots ##########
fout = open(PredictionOutputFile+".txt", 'w')
for unixtime, tau, lower, upper in zip(UnixTimes, Taus, LowerBoundaries, UpperBoundaries):
    fout.write(str(unixtime)+"\t\t")
    fout.write(str(tau)+"\t\t")
    fout.write(str(lower)+"\t\t")
    fout.write(str(upper)+"\t\t")
    fout.write("\n")
fout.close()
#######################################################################

#gOneSigmaRegion = GetRegion(UnixTimes, LowerBoundaries, UpperBoundaries)
#gOneSigmaRegion.SetFillColor(416)


#cm1 = TCanvas("cm1", "", 1500,800)
#cm1.SetMargin(0.14, 0.02, 0.14,0.02)
#cm1.cd()
#cm1.SetTickx()
#cm1.SetTicky()
#cm1.SetGridx()
#cm1.SetGridy()

#################################
## Plotting the trend if without outgassing
#################################
#TempPars = mean
#TempPars[4] = 0
#TempPars[5] = 0
#ParsForNoOutgassing, IfSth = FormPars(mean)
#pElectronLifetimeTrend.SetParameters(ParsForNoOutgassing)
#gTau_NoOutgassing = TGraph()
#for i, unixtime in enumerate(UnixTimes):
    #Tau = pElectronLifetimeTrend.GetElectronLifetime(unixtime)
    #gTau_NoOutgassing.SetPoint(i, unixtime, Tau)
#cm1.cd()
#gTau_NoOutgassing.SetLineWidth(3)
#gTau_NoOutgassing.SetLineStyle(10)
#gTau_NoOutgassing.SetLineColor(632)

######################
## Draw!
#######################


#cm1.cd()
#cm1.SetLogy()
#hframe1.Draw()
#gOneSigmaRegion.Draw("fsame")
#gEL.SetMarkerStyle(22)
#gEL.SetMarkerSize(1.5)
#gEL.SetLineWidth(3)
#gEL.SetMarkerColor(1)
#gEL.SetLineColor(1)
#gEL.Draw("psame")
#gTau.SetLineWidth(3)
#gTau.SetLineStyle(10)
#gTau.SetLineColor(1)
#gTau.Draw("lsame")
#gTau_NoOutgassing.Draw("lsame")


#cm1.cd()
#UnixTimeForCirculationIncreasing = 1462834800 + 26.*3600.*24.
#UnixTimeForGetterResume = 1466517600
#MyLine = TLine(1463.05e6, 6.96, 1463.4e6, 6.96)
#MyLine.SetLineWidth(4)
#MyLine.SetLineStyle(10)
#MyLine.SetLineColor(600)
#MyLine.DrawLine(UnixTimeForCirculationIncreasing, MinimumELifePlot, UnixTimeForCirculationIncreasing, MaximumELifePlot)
#MyLine.DrawLine(UnixTimeForGetterResume, MinimumELifePlot, UnixTimeForGetterResume, MaximumELifePlot)

#leg = TLegend(0.45, 0.15, 0.95, 0.60)
#leg.SetFillColor(0)
#leg.SetTextFont(132)
#leg.AddEntry(gOneSigmaRegion, "#pm1#sigma region of the fitted trend", "f")
#leg.AddEntry(gTau, "Best fit", "l")
##leg.AddEntry(gTau_NoOutgassing, "Best fit with no outgassing (pure exponential)", "l")
#leg.AddEntry(gEL, "E-Life data (S1 correction #tau_{s1}="+str('%.0f' % S1ExponentialConstant)+"#mus)", "pl")

#leg.Draw("same")


##################################
### SAVE !
##################################

#FigureSave = FigureOutputFile+".png"
#cm1.Print(FigureSave)
#FigureSave = FigureOutputFile+".root"
#cm1.Print(FigureSave)

#fout = open(PredictionOutputFile+".txt", 'w')
#fout.write(str(MaximumELifeRange[0])+"\t\t"+str(MaximumELifeRange[1])+"\n")
#fout.write(str(ReachingUnixtimeOf500usRange[0])+"\t\t"+str(ReachingUnixtimeOf500usRange[1])+"\n")
#fout.write(str(ReachingUnixtimeOf1000usRange[0])+"\t\t"+str(ReachingUnixtimeOf1000usRange[1])+"\n")
#fout.close()

# @ 2016-07-28
# Okay. We output the other things other than the maximum
#fout = open(PredictionOutputFile+".txt", 'w')
#for unixtime, tau, lower, upper in zip(UnixTimes, Taus, LowerBoundaries, UpperBoundaries):
#    fout.write(str(unixtime)+"\t\t")
#    fout.write(str(tau)+"\t\t")
#    fout.write(str(lower)+"\t\t")
#    fout.write(str(upper)+"\t\t")
#    fout.write("\n")
#fout.close()

#cm1.Update()
#UserInput = input("Any key: ")



