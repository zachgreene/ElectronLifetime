import glob
import os, sys
import numpy as np
import datetime, time
from Tools import *
import MongoDB_Query



if len(sys.argv) < 3:
    print('============ Usage ============')
    print('python CreateFilelists.py')
    print('<filelist output folder>')
    print('<start time>')
    print('<end time>')
    exit()

sPathToFilelists = os.path.abspath(sys.argv[1])
StartTimeStamp = sys.argv[2]
EndTimeStamp = sys.argv[3]

# unixtime from time stamps
StartUnixtime = GetUnixtime(StartTimeStamp)
EndUnixtime = GetUnixtime(EndTimeStamp)

# grab all filelists
aFilelists = glob.glob(sPathToFilelists + '/*.txt')

# access runs db
print('Accessing Runs DB')
#RunIDs, RunNames, RunSources, RunStartUnixtimes, RunEndUnixtimes = MongoDB_Query.GetRunIDsFilenamesAndSourcesFromMongo(StartTimeStamp, EndTimeStamp)
RunInfo = MongoDB_Query.GetRunIDsFilenamesAndSourcesFromMongo(StartTimeStamp, EndTimeStamp)

print('Modifying files')
# open filelists, match to runs db, and rewrite
for sFilelist in sorted(aFilelists):
    FilelistName = sFilelist.split('/')[-1]
    try:
        fin = open(sFilelist, 'r')
        lines = fin.readlines()
        fin.close()
    
        RunNamesOld, RunIDsOld, RunSourcesOld = [], [], []
        for line in lines:
            if len(line) < 2:
                continue
            line = line[:-1]
            contents = line.split('\t\t')
            assert len(contents) == 3
            RunNamesOld.append(contents[0])
            RunIDsOld.append(contents[1])
            RunSourcesOld.append(contents[2])
    except:
        print('Error: Unable to recreate %s, has %i columns' %(FilelistName, len(contents)))
        continue

    RunIDsNew = []
    RunNamesNew = []
    RunSourcesNew = []
    RunStartUnixtimesNew = []
    RunEndUnixtimesNew = []

    for RunNameOld in RunNamesOld:
        RunNamesNew.append(RunNameOld)
        RunIDsNew.append(RunInfo[RunInfo.name == RunNameOld].number.values[0])
        RunSourcesNew.append(RunInfo[RunInfo.name == RunNameOld].source.values[0])
        RunStartUnixtimesNew.append(RunInfo[RunInfo.name == RunNameOld].start.values[0])
        RunEndUnixtimesNew.append(RunInfo[RunInfo.name == RunNameOld].end.values[0])

    fout = open(sFilelist, 'w')
    for Name, ID, Source, Start, End in zip(RunNamesNew, RunIDsNew, RunSourcesNew, RunStartUnixtimesNew, RunEndUnixtimesNew):
        fout.write(Name + '\t\t' + str(int(ID)) + '\t\t' + Source + '\t\t' + str(int(Start)) + '\t\t' + str(int(End)) + '\n')
    fout.close()
    print('Successfully modified %s!' %FilelistName)
