import sys, os
from subprocess import call
import subprocess as subp
from tqdm import tqdm

if len(sys.argv) < 2:
    print('====== Usage ======')
    print('python BatchExtractLifetime.py')
    print('<list>')
    print('<lifetime output txt file>')
    print('<figure directory>')
    print('<show basic figures (0 or no, 1 for yes)>')
    print('<pax_version>')
    exit()

bShowBasicFigs = 0
pax_version = '6.6.5'

sPathToList = sys.argv[1]
LifetimeOutputFile = sys.argv[2]
sPathToFigureSaveDirectory = sys.argv[3]
if len(sys.argv) > 4:
    bShowBasicFigs = sys.argv[4]
if len(sys.argv) > 5:
    pax_version = sys.argv[5]

print('\n========    Processing electron lifetime files with pax_v' + pax_version + '    ========\n')

fin = open(sPathToList, 'r')
lines = fin.readlines()
fin.close()

if not os.path.exists(LifetimeOutputFile):
    fout = open(LifetimeOutputFile, 'w')
    fout.close() 

for line in tqdm(lines):
    sPathToFilelist = line[:-1]
    subp.call('python /home/zgreene/xenon1t/ElectronLifetime/PythonCodeV6/ExtractElectronLifetime.py ' + sPathToFilelist + ' ' + LifetimeOutputFile + ' ' + sPathToFigureSaveDirectory + ' ' + str(bShowBasicFigs) + ' ' + pax_version, shell=True)
    print('File List ' + sPathToFilelist + ' is finished\n')
