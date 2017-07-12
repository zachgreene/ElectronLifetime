import os, sys
import glob
import numpy as np
import datetime, time
from Tools import *
import MongoDB_Query

if len(sys.argv) < 3:
    print('============ Usage ============')
    print('python CreateFilelists.py')
    print('<start time>')
    print('<end time>')
    print('<duration (in hours)>')
    print('<pax_verison>')
    print('<filelist output folder (optional)>')
    print('<path to list that stores filelists (optional)>')
    exit()

StartTimeStamp = sys.argv[1]
EndTimeStamp = sys.argv[2]
Duration = float(sys.argv[3])
pax_version = sys.argv[4]

if len(sys.argv) > 5:
    sPathToOutput = sys.argv[5]
else:
    sPathToOutput = '/home/zgreene/xenon1t/ElectronLifetime/Filelists/'

if len(sys.argv) > 6:
    sPathToListOfFilelists = sys.argv[6]
else:
    sPathToListOfFilelists = '/home/zgreene/xenon1t/ElectronLifetime/Lists/List_' + EndTimeStamp[:6] + '.txt'

print('\nPath to output filelists: %s' %sPathToOutput)
print('Path to list that stores filelists: %s\n' %sPathToListOfFilelists)

##### list where minitrees are stored
MinitreePaths = ['/project/lgrandi/xenon1t/minitrees/pax_v' + pax_version + '/*Basics*']
MinitreePaths.append('/project2/lgrandi/xenon1t/minitrees/pax_v' + pax_version + '/*Basics*')
MinitreePaths.append('/project/lgrandi/zgreene/xenon1t/minitrees/pax_v' + pax_version + '/*Basics*')


#----------------------------------
# get minitrees that are on midway
#----------------------------------
ListOfMinitrees = []

for MinitreePath in MinitreePaths:
    ListOfMinitrees += glob.glob(MinitreePath)

ListOfMinitrees = ListOfMinitrees
ListOfRunNames = GetRunNamesFromMinitreeRunList(ListOfMinitrees)
ListOfRunNames = sorted(set(ListOfRunNames))



StartUnixtime = GetUnixtime(StartTimeStamp)
EndUnixtime = GetUnixtime(EndTimeStamp)

###### query db for run ids, names, and sources
print('Accessing Runs DB')
RunInfo = MongoDB_Query.GetRunIDsFilenamesAndSourcesFromMongo(StartTimeStamp, EndTimeStamp)
print(StartUnixtime, EndUnixtime)

#----------------------------------
# build filelists
#----------------------------------
print('Building filelists')

OutputFilenames = []
RunNamesToOutput, RunIDsToOutput, RunSourcesToOutput, RunStartUnixtimesToOutput, RunEndUnixtimesToOutput, RunUnixtimes = [], [], [], [], [], []


#for RunID, RunName, RunSource, RunStartUnixtime, RunEndUnixtime in sorted(zip(RunIDs, RunNames, RunSources, RunStartUnixtimes, RunEndUnixtimes), key=lambda x: x[0]): 
for i,row in RunInfo.sort_values(by='number').iterrows(): 
    RunName = row['name']
    RunID = row['number']
    RunSource = row['source']
    RunStartUnixtime = row['start']
    RunEndUnixtime = row['end']

    RunUnixtime = GetUnixtime(RunName)
    print(RunID, RunName, RunSource, end=' - ')


    if RunUnixtime < StartUnixtime:
        continue

    # make sure minitree is present
    if RunName not in ListOfRunNames:
        print('No minitree found')
        continue

    elif RunUnixtime >= StartUnixtime and len(RunNamesToOutput) == 0:
        print('List start')
        RunNamesToOutput.append(RunName)
        RunIDsToOutput.append(RunID)
        RunSourcesToOutput.append(RunSource)
        RunStartUnixtimesToOutput.append(RunStartUnixtime)
        RunEndUnixtimesToOutput.append(RunEndUnixtime)
        RunUnixtimes.append(RunUnixtime)

    elif RunUnixtime >= RunUnixtimes[0] + Duration*3600. or RunUnixtime >= EndUnixtime:
        print('list ended')
        AverageUnixtime = np.average(RunUnixtimes)        
        OutputFilename = datetime.datetime.fromtimestamp(int(AverageUnixtime)).strftime('%Y%m%d_%H%M')
        OutputFilename = OutputFilename[2:]
#        print(OutputFilename)
        AbsolutePath = os.path.abspath(sPathToOutput + '/' + OutputFilename + '_' + str(int(Duration)) + 'hrs.txt')
        print(AbsolutePath)
        fout = open(AbsolutePath, 'w') 
        for Name, ID, Source, Start, End in zip(RunNamesToOutput, RunIDsToOutput, RunSourcesToOutput, RunStartUnixtimesToOutput, RunEndUnixtimesToOutput):
            fout.write(Name + '\t\t' + str(int(ID)) + '\t\t' + Source + '\t\t' + str(int(Start)) + '\t\t' + str(int(End)) + '\n')
        fout.close()
        OutputFilenames.append(AbsolutePath)

        # if pass end time, exit loop
        if RunUnixtime >= EndUnixtime:
            print('after end')
            break

        # clear lists for next set
        RunNamesToOutput = []
        RunIDsToOutput = []
        RunSourcesToOutput = []
        RunStartUnixtimesToOutput = []
        RunEndUnixtimesToOutput = []
        RunUnixtimes = []

        RunNamesToOutput.append(RunName)
        RunIDsToOutput.append(RunID)
        RunSourcesToOutput.append(RunSource)
        RunStartUnixtimesToOutput.append(RunStartUnixtime)
        RunEndUnixtimesToOutput.append(RunEndUnixtime)
        RunUnixtimes.append(RunUnixtime)

    else:
        print('appended to list')
        RunNamesToOutput.append(RunName)
        RunIDsToOutput.append(RunID)
        RunSourcesToOutput.append(RunSource)
        RunStartUnixtimesToOutput.append(RunStartUnixtime)
        RunEndUnixtimesToOutput.append(RunEndUnixtime)
        RunUnixtimes.append(RunUnixtime)


#for RunName in ListOfRunNames:
##    print(RunName)
#    RunUnixtime = GetUnixtime(RunName)
#
#    if RunUnixtime < StartUnixtime:
#        continue
#
#    elif RunUnixtime >= StartUnixtime and len(RunNames) == 0:
#        RunNames.append(RunName)
#        RunUnixtimes.append(RunUnixtime)
#
#    elif RunUnixtime >= RunUnixtimes[0] + Duration*3600. or RunUnixtime >= EndUnixtime:
#        AverageUnixtime = np.average(RunUnixtimes)        
#        OutputFilename = datetime.datetime.fromtimestamp(int(AverageUnixtime)).strftime('%Y%m%d_%H%M')
#        OutputFilename = OutputFilename[2:]
##        print(RunNames)
##        print(OutputFilename)
#        AbsolutePath = os.path.abspath(sPathToOutput + '/' + OutputFilename + '_' + str(int(Duration)) + 'hrs.txt')
#        fout = open(AbsolutePath, 'w') 
#        for Name in RunNames:
#            fout.write(Name + '\n')
#        fout.close()
#        OutputFilenames.append(AbsolutePath)
#
#        # if pass end time, exit loop
#        if RunUnixtime >= EndUnixtime:
#            break
#
#        # clear lists for next set
#        RunNames = []
#        RunUnixtimes = []
#        RunNames.append(RunName)
#        RunUnixtimes.append(RunUnixtime)
#
#    else:
#        RunNames.append(RunName)
#        RunUnixtimes.append(RunUnixtime)

fout = open(sPathToListOfFilelists, 'w')
for OutputFilename in OutputFilenames:
   fout.write(OutputFilename + '\n')
fout.close() 

#print(OutputFilenames)
