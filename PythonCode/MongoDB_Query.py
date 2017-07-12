import os
import socket
import sys
import datetime
import time
import pandas as pd
from pymongo import MongoClient


def GetUnixTime(Filename):
    contents = Filename.split("/")
    Filename = contents[-1]
    Year = "20"+Filename[0:2]
    Month = Filename[2:4]
    Day = Filename[4:6]
    Hour = Filename[7:9]
    Minute = Filename[9:11]
    dt = datetime.datetime(int(Year), int(Month), int(Day), int(Hour), int(Minute), 0)
    return (dt, time.mktime(dt.timetuple()))



# get unixtime for beginning and end of list of runs
# accepts run names or numbers
def GetStartStopTimes(RunNamesOrNumbers):
    if isinstance(RunNamesOrNumbers, str):
        RunNamesOrNumbers = [RunNamesOrNumbers]
    RunNamesOrNumbers = sorted(RunNamesOrNumbers)

    FirstRun = RunNamesOrNumbers[0]
    LastRun = RunNamesOrNumbers[-1]

    # try to convert to number, otherwise is name
    try:
        FirstRun = int(FirstRun)
        LastRun = int(LastRun)
    except:
        pass

    StartTimeUnix, EndTimeUnix = None, None

    db_user = os.environ.get('DBUSER')
    db_pass = os.environ.get('DBPASS')

    client = MongoClient("mongodb://"+db_user+":"+db_pass+"@xenon1t-daq.lngs.infn.it:27017/run")
    collection = client['run']['runs_new']

#     get start time of first file
    query = {'$or': [{'number': FirstRun}, {'name': FirstRun}]}
    cursor = collection.find(query)
    for doc in cursor:
        StartTimeUnix = doc['trigger']['start_timestamp']

#     get end time of last file
    query = {'$or': [{'number': LastRun}, {'name': LastRun}]}
    cursor = collection.find(query)
    for doc in cursor:
        EndTimeUnix = doc['trigger']['end_trigger_processing_timestamp']

    return (StartTimeUnix, EndTimeUnix)



def GetRunIDsFilenamesAndSourcesFromMongo(StartTimeStamp, EndTimeStamp, exclude=[]):

    StartDatetime, StartUnixTime = GetUnixTime(StartTimeStamp)
    EndDatetime, EndUnixTime = GetUnixTime(EndTimeStamp)

    # Username and password for runs DB                                                    
    db_user = os.environ.get('DBUSER')
    db_pass = os.environ.get('DBPASS')

    if db_user == None or db_pass == None:
        print("You have to define a username and password for the database")
        exit()

    # Connect to database                                                                  
    client = MongoClient("mongodb://"+db_user+":"+
                         db_pass+"@xenon1t-daq.lngs.infn.it:27017/run")
    collection = client['run']['runs_new']

    # This is the query                                                                    
    query =    {
        "detector":    "tpc",
        "data": { "$elemMatch": {
#            "host": site,
#            "type": data_type,
#            "status": "transferred"
        }},
        "start": {"$gt": StartDatetime, "$lt": EndDatetime},
        "source": { "$in": [ {"type": "none"}, {"type": "Cs137"}, {"type": "Kr83m"}, {"type":"Rn220"}, {"type":"AmBe"}, {"type":"Th228"}]}
        }

    # Get the file paths and verify they exist                                             
    cursor = collection.find(query)
    Filenames = []
    RunIDs = []
    Sources = []
    RunStartUnixtimes = []
    RunEndUnixtimes = []
    for doc in cursor:
#        print(doc['name'], doc['start'])

        # remove runs that are tagged as 'bad' or 'messy'
        try:
            tags = [tag['name'] for tag in doc['tags']]
            tags = '.'.join(tags)
            if ('bad' in tags) or ('messy' in tags) or ('quake' in tags) or ('SourceMoving' in tags)
                or ('trip' in tags) or ('test' in tags) or ('crash' in tags) or ('ramping' in tags)
                or ('exception' in tags) or ('special' in tags) or ('spike' in tags):
                print(doc['name'] + ' has tags keyword - ',end='')
                print(tags)
                continue
        except:
            pass

        filename = doc['name']
        runid = doc['number']
        try:
            starttime = doc['trigger']['start_timestamp']
            endtime = doc['trigger']['end_trigger_processing_timestamp']
        except:
            print('Unable to load %s (Run number: %i)' %(filename, runid))
            continue
        UnixDatetime, unixtime = GetUnixTime(filename)
        #print(unixtime)
        if unixtime<StartUnixTime or unixtime>EndUnixTime:
            continue
        contents = filename.split("/")
        filename = contents[len(contents)-1].split(".")[0]
        Filenames.append(filename)
        RunIDs.append(runid)
        Sources.append(doc['source']['type'])
        RunStartUnixtimes.append(starttime)
        RunEndUnixtimes.append(endtime)

    dRunInfo = dict(
                    name = Filenames,
                    number = RunIDs,
                    source = Sources,
                    start = RunStartUnixtimes,
                    end = RunEndUnixtimes
                    )
#    return (RunIDs, Filenames, Sources, RunStartUnixtimes, RunEndUnixtimes)
    RunInfo =  pd.DataFrame.from_dict(dRunInfo)
    RunInfo = RunInfo[['name', 'number', 'source', 'start', 'end']]
    return RunInfo
