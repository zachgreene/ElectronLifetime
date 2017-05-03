import HistorianData
from HistorianData import *

#import ROOT
import re
import numpy as np
import scipy as sp
import os, sys
from scipy.interpolate import interp1d
import time

class MyImpurityTrend:

    def __init__(self, HistorianFile, MinUnixTime, MaxUnixTime, pars):
        Fields = [0.0583245443043, 0.0787113537639, 0.117382423495, 0.162989973122, 0.253660532825, 0.34227580980800004, 0.436223335797, 0.588573683763, 0.7941034911479999, 1.0868133019500001, 1.38497072037, 1.81605663549, 2.41539631008, 3.44978922077, 3.9223289353700004, 5.8466333187299995, 7.24083005596]
        AttachingRates = [173839448697.0, 166490987263.0, 154902864208.0, 144105038741.0, 124672767516.0, 112649545428.0, 103268822098.0, 90632673738.3, 78393350018.2, 67808393527.7, 59505871959.5, 52222275531.1, 44516323097.2, 35805292980.3, 33787097454.0, 26016257301.4, 22829772653.1]
        self.AttachingRateAsField = interp1d(Fields, AttachingRates, fill_value = 'extrapolate')
        self.ReferenceField = 0.15 # 15kV case
        self.HistorianData = MyHistorianData(HistorianFile)
        self.SetNuisanceParameters()
        self.SetTimeWindow(MinUnixTime, MaxUnixTime)
        self.SetDefaultParameters()
        self.GetGlobalHistorianData()
        self.SetParameters(pars)
        return

    def __eq__(self, other):
        # SC handler
        self.HistorianData = other.HistorianData
        self.GlobalDays = other.GlobalDays
        self.GlobalLiquidFlows = other.GlobalLiquidFlows
        self.GlobalGasFlows = other.GlobalGasFlows
        self.GlobalCoolingPowers = other.GlobalCoolingPowers
        self.GlobalCathodeVoltages = other.GlobalCathodeVoltages
        # time window
        self.MinUnixTime = other.MinUnixTime
        self.MaxUnixTime = other.MaxUnixTime
        # nuisance
        self.ReferenceField = other.ReferenceField
        self.MassGXe = other.MassGXe
        self.MassLXe = other.MassLXe
        self.LatentHeatXenon = other.LatentHeatXenon
        self.StandardXenonDensity = other.StandardXenonDensity
        self.StableCoolingPower = other.StableCoolingPower
        self.DefaultTimeStep = other.DefaultTimeStep
        # parameters
        self.InitialConcentrationGXe = other.InitialConcentrationGXe
        self.InitialConcentrationLXe = other.InitialConcentrationLXe
        self.ImpurityAttachingProbVaporization = other.ImpurityAttachingProbVaporization
        self.ImpurityAttachingProbCondensation = other.ImpurityAttachingProbCondensation
        self.OutgassingRateGXe = other.OutgassingRateGXe
        self.OutgassingRateLXe = other.OutgassingRateLXe
        self.ImpurityChangingUnixTimes = other.ImpurityChangingUnixTimes
        self.ImpurityChangingTypes = other.ImpurityChangingTypes
        self.ImpurityConcentrationChanges = other.ImpurityConcentrationChanges
        self.OutgassingRateGXeChanges = other.OutgassingRateGXeChanges # [start/end changing unixtime, changing amount], assuming linear
        self.OutgassingRateGXeDecreasingLinearConst = other.OutgassingRateGXeDecreasingLinearConst
        self.OutgassingRateLXeDecreasingLinearConst = other.OutgassingRateLXeDecreasingLinearConst
        self.OutgassingRateLXeDecreasingLinearAdditionalConsts = other.OutgassingRateLXeDecreasingLinearAdditionalConsts
        # [ start unixtime, linear coefficient]
        # interp1d
        self.inter_ConcentrationsGXe = other.inter_ConcentrationsGXe
        self.inter_ConcentrationsLXe = other.inter_ConcentrationsLXe
        self.AttachingRateAsField = other.AttachingRateAsField
        return

    def SetTimeWindow(self, MinUnixTime, MaxUnixTime):
        self.MinUnixTime = MinUnixTime
        self.MaxUnixTime = MaxUnixTime
        return

    def GetGlobalHistorianData(self):
        Npoints = (self.MaxUnixTime - self.MinUnixTime) / 3600. / 24. / self.DefaultTimeStep + 1
        self.GlobalDays = np.linspace(
                                         (self.MinUnixTime - self.HistorianData.GetReferenceUnixTime())/3600./24.,
                                         (self.MaxUnixTime - self.HistorianData.GetReferenceUnixTime())/3600./24.,
                                         int(Npoints)
                                        )
        self.GlobalLiquidFlows = []
        self.GlobalGasFlows = []
        self.GlobalCoolingPowers = []
        self.GlobalCathodeVoltages = []
        for day in self.GlobalDays:
            unixtime = day*3600.*24. + self.HistorianData.GetReferenceUnixTime()
            LiquidFlow, GasFlow, CoolingPower, CathodeVoltage = self.HistorianData.GetHistorian(unixtime)
            self.GlobalLiquidFlows.append(LiquidFlow)
            self.GlobalGasFlows.append(GasFlow)
            self.GlobalCoolingPowers.append(CoolingPower)
            self.GlobalCathodeVoltages.append(CathodeVoltage)
        return

    # update the global historian data because the getter deficiency
    # It must be after running self.GetGlobalHistorianData()
    def UpdateGlobalHistorianData(self):
        # to save time
        # only update the period when there assumes to be a getter deficiency
        Npoints = (self.MaxUnixTime - self.MinUnixTime) / 3600. / 24. / self.DefaultTimeStep + 1
        Configs = self.HistorianData.GetGetterDeficiencyConfigs()
        for Config in Configs:
            StartUnixtime = Config[0]
            EndUnixtime = Config[1]
            StartIndex = int((StartUnixtime - self.MinUnixTime) / 3600. / 24. / self.DefaultTimeStep) + 1
            EndIndex = int((EndUnixtime - self.MinUnixTime) / 3600. / 24. / self.DefaultTimeStep+1)
            for index in range(StartIndex, EndIndex):
                if index>=Npoints:
                    continue
                unixtime = self.GlobalDays[index]*3600.*24.+self.HistorianData.GetReferenceUnixTime()
                LiquidFlow, GasFlow, CoolingPower, CathodeVoltage = self.HistorianData.GetHistorian(unixtime)
                self.GlobalLiquidFlows[index] = LiquidFlow
                self.GlobalGasFlows[index] = GasFlow
        return

    def SetNuisanceParameters(self):
        # totally 4 nuisance parameters are needed
        # the total mass in the three relevant region
        # and the latent heat of Xenon
        # plus the necessary stand xenon density
        # pluse the default integration step (1hr)
        self.MassGXe = 23. # kg
        self.MassLXe = 401.+(2756.- 0.) # kg. The LXe mass total
        self.LatentHeatXenon = 95.587e3 # J/kg
        self.StandardXenonDensity = 5.894 # g/L
        self.StableCoolingPower = 140. #W
        self.DefaultTimeStep = 1./24. # days (1hr)
        return

    def SetDefaultParameters(self):
        self.InitialConcentrationGXe = 100. # ppb
        self.InitialConcentrationLXe = 100. # ppb
        self.ImpurityAttachingProbVaporization = 1. # 100%
        self.ImpurityAttachingProbCondensation = 1. # 100%
        self.OutgassingRateGXe = 1. # kg*ppb/day
        self.OutgassingRateLXe = 1. # kg*ppb/day
        self.ImpurityChangingUnixTimes = [] #
        self.ImpurityChangingTypes = [] # type 0 in gas, 1 in liquid, 2  in the backing purified gas
        self.ImpurityConcentrationChanges = [] # mol
        self.OutgassingRateGXeChanges = []
        self.OutgassingRateGXeDecreasingLinearConst = 10000. # days
        self.OutgassingRateLXeDecreasingLinearConst = 10000. # days
        self.OutgassingRateLXeDecreasingLinearAdditionalConsts = []
        return

    def SetParameters(self, pars):
        # there're 12 parameters needed
        # also will start the calculation if the parameters have been changed
        # @2016-12-12 either 12 or 13 pars
        # @2017-02-07 either 13 or 14 pars
        # @2017-02-13 either 14 or 15 pars
        if not (len(pars)==15 or len(pars)==16):
            raise ValueError("Number of parameters not enough!")
        IfSameAsPrevious = self.CheckIfSame(pars)
        self.InitialConcentrationGXe = pars[0]
        self.InitialConcentrationLXe = pars[1]
        if self.InitialConcentrationGXe<0:
            self.InitialConcentrationGXe = 0.
        if self.InitialConcentrationLXe<0:
            self.InitialConcentrationLXe = 0.
        self.ImpurityAttachingProbVaporization = pars[2]
        if self.ImpurityAttachingProbVaporization<0:
            self.ImpurityAttachingProbVaporization=0
        if self.ImpurityAttachingProbVaporization>1:
            self.ImpurityAttachingProbVaporization = 1.
        self.ImpurityAttachingProbCondensation = pars[3]
        if self.ImpurityAttachingProbCondensation<0:
            self.ImpurityAttachingProbCondensation=0
        if self.ImpurityAttachingProbCondensation>1:
            self.ImpurityAttachingProbCondensation = 1.
        self.OutgassingRateGXe = pars[4]
        if self.OutgassingRateGXe<0:
            self.OutgassingRateGXe=0
        self.OutgassingRateLXe = pars[5]
        if self.OutgassingRateLXe<0:
            self.OutgassingRateLXe=0
        self.ImpurityChangingUnixTimes = pars[6]
        self.ImpurityChangingTypes = pars[7]
        self.ImpurityConcentrationChanges = pars[8]
        self.OutgassingRateGXeChanges = pars[9]
        self.OutgassingRateLXeDecreasingLinearConst = pars[11]
        self.OutgassingRateGXeDecreasingLinearConst = pars[10]
        if self.OutgassingRateGXeDecreasingLinearConst<0.001/3600./24.:
            # cannot be less than 1ms
            self.OutgassingRateGXeDecreasingLinearConst = 0.001/3600./24.
        if self.OutgassingRateLXeDecreasingLinearConst<0.001/3600./24.:
            # cannot be less than 1ms
            self.OutgassingRateLXeDecreasingLinearConst = 0.001/3600./24.
        self.OutgassingRateLXeDecreasingLinearAdditionalConsts = pars[12]
        if len(pars)==14:
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[13])
        if len(pars)==15:
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[13])
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[14])
        if len(pars)==16:
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.PopOneGetterDeficiencyConfig()
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[13])
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[14])
            self.HistorianData.AddOneGetterDeficiencyConfig(pars[15])
        if not IfSameAsPrevious:
            self.CalculateImpurityConcentration()
        return

    def GetParameters(self):
        return [
                                    self.InitialConcentrationGXe,
                                    self.InitialConcentrationLXe,
                                    self.ImpurityAttachingProbVaporization,
                                    self.ImpurityAttachingProbCondensation,
                                    self.OutgassingRateGXe,
                                    self.OutgassingRateLXe,
                                    self.ImpurityChangingUnixTimes,
                                    self.ImpurityChangingTypes,
                                    self.ImpurityConcentrationChanges,
                                    self.OutgassingRateGXeChanges,
                                    self.OutgassingRateLXeDecreasingLinearConst,
                                    self.OutgassingRateGXeDecreasingLinearConst,
                                    self.OutgassingRateLXeDecreasingLinearAdditionalConsts,
                                    self.HistorianData.GetGetterDeficiencyConfigs(),
                                   ]

    def CheckIfSame(self, pars):
        PreviousPars = self.GetParameters()
        for previous_par, par in zip(PreviousPars, pars):
            if not previous_par==par:
                return False
        return

    def GetHistorianData(self):
        return self.HistorianData

    def CalculateImpurityConcentration(self):
        # main function for calculating the impurity trend
        # print(str(self.MinUnixTime) + " <-> " + str(self.MaxUnixTime) )
        Npoints = len(self.GlobalDays)
        ConcentrationsGXe = []
        ConcentrationsLXe = []
        TrueTimeStep = (self.MaxUnixTime - self.MinUnixTime) / 3600. / 24. / float(Npoints - 1)
        self.UpdateGlobalHistorianData()
        for i, (day, LiquidFlow, GasFlow, CoolingPower, CathodeVoltage) in enumerate(zip(
                self.GlobalDays, 
                self.GlobalLiquidFlows,
                self.GlobalGasFlows,
                self.GlobalCoolingPowers,
                self.GlobalCathodeVoltages,
                )):
            if i==0:
                ConcentrationsGXe.append(self.InitialConcentrationGXe)
                ConcentrationsLXe.append(self.InitialConcentrationLXe)
                continue
            unixtime = day*3600.*24. + self.HistorianData.GetReferenceUnixTime()
            # Previous concentrations
            PreviousConcentrationGXe = ConcentrationsGXe[i-1]
            PreviousConcentrationLXe = ConcentrationsLXe[i-1]
            # differential concentration of GXe & LXe
            ConcentrationChangeGXe = self.OutgassingRateGXe * (1. -  day / self.OutgassingRateGXeDecreasingLinearConst )# linear
            # define LXe decreasing rate
            theOutgassingRateLXeDecreasingLinearConst = self.OutgassingRateLXeDecreasingLinearConst 
            for i, config in enumerate(self.OutgassingRateLXeDecreasingLinearAdditionalConsts):
                if unixtime < config[0]:
                    break
                theOutgassingRateLXeDecreasingLinearConst = config[1]
            ConcentrationChangeLXe = self.OutgassingRateLXe * (1. - day / theOutgassingRateLXeDecreasingLinearConst ) # linear
            # check if necessary to change the outgassing level in GXe
            for changes in self.OutgassingRateGXeChanges:
                if len(changes)<3:
                    raise ValueError("Change array number < 3")
                    continue
                unixtime_start = changes[0]
                unixtime_end = changes[1]
                change_amount = changes[2]
                if unixtime_start>=unixtime_end:
                    continue
                if unixtime<unixtime_start:
                    continue
                elif unixtime<unixtime_end:
                    # using exponential model
                    starting_outgassing = self.OutgassingRateGXe * (1. - (unixtime_start - self.HistorianData.GetReferenceUnixTime()) / 24. / 3600. / self.OutgassingRateGXeDecreasingLinearConst )
                    if starting_outgassing<=0:
                        break
                    if starting_outgassing+change_amount <= 0:
                        ConcentrationChangeGXe = 0
                        break
                    change_exp = (unixtime_end - unixtime_start) / (-np.log( (starting_outgassing + change_amount) / starting_outgassing )) # in sec
                    ConcentrationChangeGXe -= starting_outgassing * ( 1. - np.exp( - (unixtime - unixtime_start) / change_exp ) )
                else:
                    ConcentrationChangeGXe += change_amount
            if ConcentrationChangeGXe<0:
                #print("ConcentrationGXe < 0 and set to 0")
                ConcentrationChangeGXe = 0
            ###################
            # the individual terms
            # all in positive 
            ##################
            # Gas flow term
            GasFlowTerm = GasFlow*self.StandardXenonDensity*PreviousConcentrationGXe
            GasFlowTerm *= 1e-3*60.*24. # to kg*ppb / day
            # Liquid flow term
            LiquidFlowTerm = LiquidFlow*self.StandardXenonDensity*PreviousConcentrationLXe
            LiquidFlowTerm *= 1e-3*60.*24.
            # Vaporization
            VaporizationTerm = self.ImpurityAttachingProbVaporization*self.StableCoolingPower*PreviousConcentrationLXe/self.LatentHeatXenon
            VaporizationTerm *= 3600.*24. # to kg*ppb/day
            # Condenzation
            CondensationTerm = self.ImpurityAttachingProbCondensation*CoolingPower*PreviousConcentrationGXe/self.LatentHeatXenon
            CondensationTerm *= 3600.*24. # to kg*ppb/day
            # concentratin change
            ConcentrationChangeGXe += -GasFlowTerm
            ConcentrationChangeGXe += VaporizationTerm
            ConcentrationChangeGXe += -CondensationTerm
            ConcentrationChangeLXe += -LiquidFlowTerm
            ConcentrationChangeLXe += -VaporizationTerm
            ConcentrationChangeLXe += CondensationTerm
            # time step & mass normalization
            ConcentrationChangeGXe *= TrueTimeStep/self.MassGXe
            ConcentrationChangeLXe *= TrueTimeStep/self.MassLXe
            # check if need the change
            for i, unixtime_change in enumerate(self.ImpurityChangingUnixTimes):
                if unixtime>unixtime_change and unixtime-TrueTimeStep*3600.*24. < unixtime_change:
                    if self.ImpurityChangingTypes[i]==0:
                        ConcentrationChangeGXe += self.ImpurityConcentrationChanges[i] / self.MassGXe *0.133 * 1e9
                    elif self.ImpurityChangingTypes[i]==1:
                        ConcentrationChangeLXe += self.ImpurityConcentrationChanges[i] / self.MassLXe * 0.133 * 1e9
                    else:
                        FractionGoToGas = (CoolingPower - self.StableCoolingPower) / self.LatentHeatXenon * 1.e3 * 60. / ((LiquidFlow+GasFlow)*self.StandardXenonDensity)
                        if FractionGoToGas>1:
                            FractionGoToGas = 1.
                        elif FractionGoToGas<0:
                            FractionGoToGas=0.
                        ConcentrationChangeGXe += FractionGoToGas*self.ImpurityConcentrationChanges[i] / self.MassGXe * 0.133 * 1e9
                        ConcentrationChangeLXe += (1. - FractionGoToGas)*self.ImpurityConcentrationChanges[i] / self.MassLXe * 0.133 * 1e9
                    break
            # append
            ConcentrationsGXe.append(PreviousConcentrationGXe + ConcentrationChangeGXe)
            ConcentrationsLXe.append(PreviousConcentrationLXe + ConcentrationChangeLXe)
        self.inter_ConcentrationsGXe = interp1d(self.GlobalDays, ConcentrationsGXe)
        self.inter_ConcentrationsLXe = interp1d(self.GlobalDays, ConcentrationsLXe)
        return

    def GetConcentrations(self, unixtime):
        return [self.GetGXeConcentration(unixtime),
                     self.GetLXeConcentration(unixtime),
                    ]

    def GetGXeConcentration(self, unixtime):
        EffUnixTime = unixtime
        if unixtime<self.MinUnixTime:
            EffUnixTime = self.MinUnixTime
        if unixtime>self.MaxUnixTime:
            EffUnixTime = self.MaxUnixTime
        EffDay = (EffUnixTime - self.HistorianData.GetReferenceUnixTime()) / 3600. / 24.
        return self.inter_ConcentrationsGXe(EffDay)

    def GetLXeConcentration(self, unixtime):
        EffUnixTime = unixtime
        if unixtime<self.MinUnixTime:
            EffUnixTime = self.MinUnixTime
        if unixtime>self.MaxUnixTime:
            EffUnixTime = self.MaxUnixTime
        EffDay = (EffUnixTime - self.HistorianData.GetReferenceUnixTime()) / 3600. / 24.
        return self.inter_ConcentrationsLXe(EffDay)

    def GetStandardXenonDensity(self):
        return self.StandardXenonDensity

    def GetStableCoolingPower(self):
        return self.StableCoolingPower

    def GetLatentHeat(self):
        return self.LatentHeatXenon

    def GetGXeMass(self):
        return self.MassGXe

    def GetLXeMass(self):
        return self.MassLXe

    def GetDefaultTimeStep(self):
        return self.DefaultTimeStep

    # new @ 2016-10-29
    def GetAttachingRateCorrectionFactor(self, unixtime):
        _, _, _, CathodeVoltage = self.HistorianData.GetHistorian(unixtime)
        Field = CathodeVoltage / 100. 
        return self.AttachingRateAsField (Field) / self.AttachingRateAsField (self.ReferenceField)

    # new @ 2016-11-16
    def ManualChangingFlows(self, NewLiquidFlow, NewGasFlow, start_unixtime, end_unixtime):
        self.HistorianData.AddGasOnlyPeriod([start_unixtime, end_unixtime])
        for i, day in enumerate(self.GlobalDays):
            unixtime = self.HistorianData.GetReferenceUnixTime() + day * 24. * 3600.
            if unixtime<start_unixtime or unixtime>end_unixtime:
                continue
            self.GlobalLiquidFlows[i] = NewLiquidFlow
            self.GlobalGasFlows[i] = NewGasFlow
        return

