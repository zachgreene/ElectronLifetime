import ImpurityTrend
from ImpurityTrend import *

#import ROOT
import re
import numpy as np
import scipy as sp
import os, sys
from scipy.interpolate import interp1d

class MyElectronLifetimeTrend:

    def __init__(self, HistorianFile, MinUnixTime, MaxUnixTime, pars):
        # need 1 par3 for electron lifetime trend
        # need 8 pars for impurity trend
        # 2016-12-12 need either 13 or 14 parss
        if not (len(pars)==15 or len(pars)==16 or len(pars)==17):
            raise ValueError("Parameters for electron lifetime are not enough!")
        pars_electronlifetime = pars[0:1]
        pars_impurity = pars[1:]
        self.ImpurityTrend = MyImpurityTrend(HistorianFile, MinUnixTime, MaxUnixTime, pars_impurity)
        self.SetDefaultElectronLifetimeParameters(pars_electronlifetime)
        self.DefaultTimeStep = 1./24. # 1 hour
        return

    def SetDefaultElectronLifetimeParameters(self, pars):
        if not len(pars)==1:
            raise ValueError("Parameters for default electron lifetime setting not enough!")
        self.ImpurityAttachingRate = pars[0]
        return

    def SetParameters(self, pars):
        if not (len(pars)==15 or len(pars)==16 or len(pars)==17):
            raise ValueError("Parameters are not enough for electron lifetime!")
        self.ImpurityAttachingRate = np.abs(pars[0])
        if self.ImpurityAttachingRate<1e-40:
            self.ImpurityAttachingRate = 1e-40
        pars_impurity = pars[1:]
        self.ImpurityTrend.SetParameters(pars_impurity)
        return    

    def GetParameters(self):
        pars =  [
                      self.ImpurityAttachingRate
                     ]
        pars.extend(self.ImpurityTrend.GetParameters())
        return pars

    def GetImpurityTrend(self):
        return self.ImpurityTrend

    def GetElectronLifetime(self, unixtime):
        # the electron lifetime trend model
        Ig, Il = self.ImpurityTrend.GetConcentrations(unixtime)
        AttachingRateCorrection = self.ImpurityTrend.GetAttachingRateCorrectionFactor(unixtime)
        if Il==0 or self.ImpurityAttachingRate==0 or AttachingRateCorrection==0:
            return np.inf
        return 1./self.ImpurityAttachingRate/AttachingRateCorrection/Il

    def ManualChangingFlows(self, NewLiquidFlow, NewGasFlow, start_unixtime, end_unixtime):
        self.ImpurityTrend.ManualChangingFlows(NewLiquidFlow, NewGasFlow, start_unixtime, end_unixtime)
        return
