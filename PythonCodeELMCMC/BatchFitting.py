import subprocess as subp
import sys, os


if len(sys.argv)<3:
    print('======== Syntax ==========')
    print('python BatchFitting.py .......')
    print('<SC file>')
    print('<Output pickle result>')
    print('<Fit data txt>')
    print('<Rn Fit data txt>')
    print('<number of walkers>')
    print('<number of iteractions>')
    print('<SubmitID>')
    print('<(optional) 1 for public node, 0 for not (default)>')
    print('<(optional) input pickle pre-result>')
    exit()

bPublicNode = 0
InputPickleFile = ''

HistorianFile = sys.argv[1]
OutputPickleFile = sys.argv[2]
ElectronLifetimeDataFile = sys.argv[3]
RnElectronLifetimeDataFile = sys.argv[4]
NumWalkers = sys.argv[5]
NumIterations = sys.argv[6]
SubmitID = sys.argv[7]
if len(sys.argv) > 8:
    bPublicNode = int(sys.argv[8])
if len(sys.argv) > 9:
    InputPickleFile = sys.argv[9]


# get absolute paths to files
HistorianFile = os.path.abspath(HistorianFile)
OutputPickleFile = os.path.abspath(OutputPickleFile)
ElectronLifetimeDataFile = os.path.abspath(ElectronLifetimeDataFile)
RnElectronLifetimeDataFile = os.path.abspath(RnElectronLifetimeDataFile)
if len(sys.argv) > 9:
    InputPickleFile = os.path.abspath(InputPickleFile)

OutputPredictionFile = OutputPickleFile.split('/')[-1].split('_')[-1].split('.')[0]
OutputPredictionFile = '/'.join(OutputPickleFile.split('/')[:-2]) + '/TXTs/Prediction_' + OutputPredictionFile
BurnInWalkers = str(int(NumIterations) - 200)


EXE1 = '/home/zgreene/xenon1t/ElectronLifetime/PythonCodeELMCMC/FitElectronLifetime.py'
ARGS1 = ' '.join([HistorianFile, OutputPickleFile, ElectronLifetimeDataFile, RnElectronLifetimeDataFile, NumWalkers, NumIterations, InputPickleFile])

EXE2 = '/home/zgreene/xenon1t/ElectronLifetime/PythonCodeELMCMC/PredictElectronLifetime.py'
ARGS2 = ' '.join([HistorianFile, OutputPickleFile, OutputPredictionFile, BurnInWalkers])


#print('python ' + EXE1 + ' ' + ARGS1 + '\n')
#print('python ' + EXE2 + ' ' + ARGS2 + '\n')
#exit()

# confirm BatchFitting.py and FitElectronLifetime.py are from same version
#BatchFitPath = os.path.abspath(sys.argv[0])
#ContentsBatchFit = BatchFitPath.split('/')
#ContentsFit = EXE1.split('/')

#if ContentsBatchFit[-2] != ContentsFit[-2]:
#    print('\nERROR: %s/BatchFitting.py and %s/FitEelectronLifetime.py are from different versions, exiting\n' %(ContentsBatchFit[-2], ContentsFit[-2]))
#    exit()


# create submit file
SubmitPath = '/home/zgreene/xenon1t/ElectronLifetime/JobSubmit/' + SubmitID
if os.path.exists(SubmitPath):
    subp.call("rm -r " + SubmitPath, shell=True)
subp.call("mkdir " + SubmitPath, shell=True)

SubmitFile = SubmitPath + '/submit'
if os.path.exists(SubmitFile):
    subp.call("rm " + SubmitFile, shell=True)



subp.call("echo '#!/bin/bash\n' >>" + SubmitFile, shell=True)
subp.call("echo '#SBATCH --output=" + SubmitPath + "/myout_" + str(SubmitID) + ".txt' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --error=" + SubmitPath + "/myerr_" + str(SubmitID) + ".txt' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --time=24:00:00\n' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --account=pi-lgrandi\n' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --nodes=1\n' >> " + SubmitFile, shell=True)
if not bPublicNode:
    subp.call("echo '#SBATCH --qos=xenon1t\n' >> " + SubmitFile, shell=True)
    subp.call("echo '#SBATCH --partition=xenon1t\n' >> " + SubmitFile, shell=True)
    subp.call("echo '#SBATCH --ntasks-per-node=28\n' >> " + SubmitFile, shell=True)
else:
    subp.call("echo '#SBATCH --qos=xenon1t-kicp\n' >> " + SubmitFile, shell=True)
    subp.call("echo '#SBATCH --partition=kicp\n' >> " + SubmitFile, shell=True)
    subp.call("echo '#SBATCH --ntasks-per-node=16\n' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --mail-type=END\n' >> " + SubmitFile, shell=True)
subp.call("echo '#SBATCH --mail-user=zgreene\n' >> " + SubmitFile, shell=True)
subp.call("echo '. /home/zgreene/ENV/GlobalPAXEnv.sh\n' >> " + SubmitFile, shell=True)
#subp.call("echo 'echo $PYTHONPATH\n' >> " + SubmitFile, shell=True)
#subp.call("echo 'echo $PATH\n' >> " + SubmitFile, shell=True)
#subp.call("echo '. /home/zgreene/ENV/BatchReductionEnv.sh\n' >> " + SubmitFile, shell=True)
subp.call("echo 'python " + EXE1 + " " + ARGS1 +"\n' >> " + SubmitFile, shell=True)
subp.call("echo 'python " + EXE2 + " " + ARGS2 +"\n' >> " + SubmitFile, shell=True)


#subp.call("python " + EXE + " " + ARGS, shell=True)

# submit batch
subp.call("cd " + SubmitPath + ";sbatch " + SubmitFile + "; cd -", shell=True)
