import numpy as np
import pandas as pd
import datetime
import sys

import FormPars
import MCMC_Tools

if len(sys.argv) < 3:
    print("======== Syntax: =======")
    print("python ConvertEvolution.py .....")
    print("<Rn222 elife data txt file>")
#    print("<Kr83 elife data txt file>")
    print("<prediction txt file> ")
    exit()


PathToLifetimesRn222      = sys.argv[1]
#PathToLifetimesKr83       = sys.argv[2]
PathToLifetimePredictions = sys.argv[2]


PathContents = PathToLifetimePredictions.split('/')
PathToKrLifetimeFile = '/'.join(PathContents[:-1]) + '/Kr83m/Kr_' + PathContents[-1]

ChangeVal, ChangeValErr = FormPars.GetKrCorrection()

(
    Rn222Unixtimes,
    Rn222UnixtimeErrors,
    Rn222ELifeValues,
    Rn222ELifeValueErrors) = np.asarray(MCMC_Tools.LoadFitData('Rn222', PathToFile=PathToLifetimesRn222))

#(
#    KrUnixtimes,
#    KrUnixtimeErrors,
#    KrELifeValues,
#    KrELifeValueErrors) = np.asarray(MCMC_Tools.LoadFitData('Kr83', PathToFile=PathToLifetimesKr83))

ScienceRunStartUnixtime = 1486054320
ScienceRunEndUnixtime   = max(Rn222Unixtimes)
ScienceRunStartDatetime = datetime.datetime.fromtimestamp(ScienceRunStartUnixtime)
ScienceRunEndDatetime   = datetime.datetime.fromtimestamp(ScienceRunEndUnixtime)

LastPointUnixtime       = ScienceRunEndUnixtime
DaysAfterLastPoint      = 100

# load Rn222 prediction
(
    PredictionUnixtimes,
    Rn222PredictedELifes,
    Rn222PredictedELifeLows,
    Rn222PredictedELifeUps,
    Rn222PredictedELifeLowErrs,
    Rn222PredictedELifeUpErrs) = MCMC_Tools.LoadPredictions(PathToLifetimePredictions,
                                                            LastPointUnixtime=LastPointUnixtime,
                                                            DaysAfterLastPoint=DaysAfterLastPoint)


# change to Kr
(
    KrPredictedELifes,
    KrPredictedELifeLows,
    KrPredictedELifeUps,
    KrPredictedELifeLowErrs,
    KrPredictedELifeUpErrs) = MCMC_Tools.ChangeElectronLifetime(Rn222PredictedELifes,
                                                                Rn222PredictedELifeLows,
                                                                Rn222PredictedELifeUps,
                                                                Rn222PredictedELifeLowErrs,
                                                                Rn222PredictedELifeUpErrs,
                                                                ChangeVal,
                                                                ChangeValErr
                                                                )

fout = open(PathToKrLifetimeFile, 'w')
for unixtime, tau, lower, upper in zip(PredictionUnixtimes,
                                        KrPredictedELifes,
                                        KrPredictedELifeLows,
                                        KrPredictedELifeUps
                                        ):
    fout.write(str(unixtime)+"\t\t")
    fout.write(str(tau)+"\t\t")
    fout.write(str(lower)+"\t\t")
    fout.write(str(upper)+"\t\t")
    fout.write("\n")
fout.close()
