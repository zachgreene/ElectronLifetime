import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
#from Cuts import EvalCut
import hax


def GetFilelistNameFromDirectory(sPathToFigureSaveDirectory):
    if sPathToFigureSaveDirectory[-1] == '/':
        sPathToFigureSaveDirectory = sPathToFigureSaveDirectory[:-1] 
    FilelistName = sPathToFigureSaveDirectory.split('/')[-1]
    return FilelistName


def PlotS1S1Asym(data, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))
    plt.scatter(data.s1, data.s1_aft,
                c=data.z, marker='.', edgecolors=None,
                vmin=-100, vmax=0, linewidths=0, s=6, cmap='viridis')
    plt.xlim(0,100000)
    plt.ylim(0,0.4)
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S1 Asymmetry (no cuts)', fontsize=20)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S1 Area Fraction Top')
    plt.colorbar(label='z [cm]')

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S1S1AsymNoCuts_' + FilelistName + '.png')
    return fig


def PlotS1S2(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))
    S1Bins = np.linspace(0,1.2*CutLimits['UpLimitS1'],200)
    S2Bins = np.linspace(0,1.2*CutLimits['UpLimitS2'],200)
    plt.hist2d(data.s1, data.s2,
           bins=[S1Bins,S2Bins],
           norm=colors.LogNorm(),cmap='viridis')

    # box in region of acceptance
    plt.hlines([CutLimits['LowLimitS2'], CutLimits['UpLimitS2']],
               CutLimits['LowLimitS1'], CutLimits['UpLimitS1'], 
               color='r', linestyle='solid')
    plt.vlines([CutLimits['LowLimitS1'], CutLimits['UpLimitS1']],
               CutLimits['LowLimitS2'], CutLimits['UpLimitS2'],
               color='r', linestyle='solid')

    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S2 (no cuts)', fontsize=20)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S2 [PE]')

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S1S2NoCuts_' + FilelistName + '.png')
    return fig



def PlotS1S2WithCuts(data, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))

    # draw events that fail cuts
#    plt.scatter(data[Xcuts==0]['s1'].values, data[Xcuts==0]['s2'].values,
#                c=data['z'], marker='.', edgecolors=None,
#                vmin=-100, vmax=0, alpha=0.5, linewidths=0, s=6, cmap='viridis')

    # draw events that meet cuts
    plt.scatter(data.s1, data.s2,
                c=data.z, cmap='viridis', vmin=-100, vmax=0,
                marker='.', edgecolors=None, linewidths=0, s=12)

    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S2 (all cuts)', fontsize=20)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S2 [PE]')
    plt.colorbar(label='z [cm]')

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S1S2AllCuts_' + FilelistName + '.png')
    return fig


def PlotS1Asym(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))

    # draw scatter plot without cuts first
    plt.subplot(2, 1, 1)
    plt.scatter(data.s1, data.s2, c=data.s1_aft,
                     marker='.', edgecolors=None, vmin=0.0, vmax=0.4, linewidths=0, s=4,cmap='viridis')

    plt.xlim(20000,80000)
    plt.ylim(0,150000)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S2 [PE]')
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S2 (no cuts)', fontsize=20)
    plt.colorbar(label='S1 Asymmetry')

    # draw histograms with and without cuts
    plt.subplot(2, 1, 2)
    plt.hist(data.s1_aft, bins=200, range=[0,1], log=True, label='No Cuts', edgecolor='none')

    # get cuts
    data = hax.cuts.range_selection(data, 's1', (CutLimits['LowLimitS1'], CutLimits['UpLimitS1']))
    data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))

    plt.hist(data.s1_aft, bins=200, range=[0,1], log=True, color='r', label='S1 & S2 Cuts',
            edgecolor='none')

    plt.vlines([CutLimits['LowLimitAsymS1'], CutLimits['UpLimitAsymS1']],
               plt.ylim()[0], plt.ylim()[1], color='mediumturquoise', linewidth=2, label='S1 Asymmetry cuts')

    plt.annotate(s='', xy=(CutLimits['LowLimitAsymS1'],0.7*plt.ylim()[1]),
                 xytext=(CutLimits['UpLimitAsymS1'],0.7*plt.ylim()[1]), size=20,
                 arrowprops=dict(arrowstyle='<|-|>', shrinkA=2, shrinkB=2, fc='mediumturquoise', ec='mediumturquoise'))

    plt.xlabel('S1 Asymmetry')
    plt.ylabel('Counts')
    plt.legend()

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S1Asym_' + FilelistName + '.png')
    return fig


def PlotS2Asym(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))

    # draw scatter plot without cuts first
    plt.subplot(2, 1, 1)
    plt.scatter(data.s1, data.s2, c=data.s2_aft,
                     marker='.', edgecolors=None, vmin=0.3, vmax=0.75, linewidths=0, s=4,cmap='viridis')

    plt.xlim(20000,80000)
    plt.ylim(0,150000)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S2 [PE]')
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S2 (no cuts)', fontsize=20)
    plt.colorbar(label='S2 Asymmetry')

    # draw histograms with and without cuts
    plt.subplot(2, 1, 2)
    plt.hist(data.s2_aft, bins=200, range=[0,1], log=True, label='No Cuts', edgecolor='none')

    # get cuts
    data = hax.cuts.range_selection(data, 's1', (CutLimits['LowLimitS1'], CutLimits['UpLimitS1']))
    data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))

    plt.hist(data.s2_aft, bins=200, range=[0,1], log=True, color='r', label='S1 & S2 Cuts',
            edgecolor='none')

    plt.vlines([CutLimits['LowLimitAsymS2'], CutLimits['UpLimitAsymS2']],
               plt.ylim()[0], plt.ylim()[1], color='mediumturquoise', linewidth=2, label='S2 Asymmetry cuts')

    plt.annotate(s='', xy=(CutLimits['LowLimitAsymS2'],0.7*plt.ylim()[1]),
                 xytext=(CutLimits['UpLimitAsymS2'],0.7*plt.ylim()[1]), size=20,
                 arrowprops=dict(arrowstyle='<|-|>', shrinkA=2, shrinkB=2, fc='mediumturquoise', ec='mediumturquoise'))

    plt.xlabel('S2 Asymmetry')
    plt.ylabel('Counts')
    plt.legend()

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S2Asym_' + FilelistName + '.png')
    return fig


#def PlotS2Asym(data, cuts, sPathToFigureSaveDirectory=''):
#    fig = plt.figure(figsize=(10,8))
#
#    # draw scatter plot without cuts first
#    plt.subplot(2, 1, 1)
#    plt.scatter(data['s1'], data['s2'], c=data['s2_aft'],
#                     marker='.', edgecolors=None, vmin=0.3, vmax=0.75, linewidths=0, s=4,cmap='viridis')
#
#    plt.xlim(20000,80000)
#    plt.ylim(0,150000)
#    plt.xlabel('S1 [PE]')
#    plt.ylabel('S2 [PE]')
#    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
#    plt.title(FilelistName + ': S1 vs. S2 (no cuts)', fontsize=20)
#    plt.colorbar(label='S2 Asymmetry')
#
#    # get cuts
#    Xs1 = EvalCut(data, cuts, 'Xs1')
#    Xs2 = EvalCut(data, cuts, 'Xs2')
#
#    # draw histograms with and without cuts
#    plt.subplot(2, 1, 2)
#    plt.hist(data['s2_aft'], bins=200, range=[0,1], log=True, label='No Cuts', edgecolor='none')
#    plt.hist(data[Xs1 & Xs2]['s2_aft'], bins=200, range=[0,1], log=True, color='r', label='S1 & S2 Cuts',
#            edgecolor='none')
#
#    plt.vlines([cuts['LowLimitAsymS2'], cuts['UpLimitAsymS2']],
#               plt.ylim()[0], plt.ylim()[1], color='mediumturquoise', linewidth=2, label='S1 Asymmetry cuts')
#
#    plt.annotate(s='', xy=(cuts['LowLimitAsymS2'],0.7*plt.ylim()[1]),
#                 xytext=(cuts['UpLimitAsymS2'],0.7*plt.ylim()[1]), size=20,
#                 arrowprops=dict(arrowstyle='<|-|>', shrinkA=2, shrinkB=2, fc='mediumturquoise', ec='mediumturquoise'))
#
#    plt.xlabel('S2 Asymmetry')
#    plt.ylabel('Counts')
#    plt.legend()
#
#    plt.tight_layout()
#    if len(sPathToFigureSaveDirectory) > 1:
#        plt.savefig(sPathToFigureSaveDirectory + 'S2Asym_' + FilelistName + '.png')
#    return fig
#

def PlotRadius(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    data = hax.cuts.range_selection(data, 's1', (CutLimits['LowLimitS1'], CutLimits['UpLimitS1']))
    data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))

    plt.scatter(data.x, data.y,
                c=data.dt,
                marker='.', linewidths=0, s=6, vmin=5, vmax=CutLimits['MaxValueDt'], cmap='viridis')
    plt.xlim(-50,50)
    plt.ylim(-50,50)
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': X vs. Y (S1 & S2 cuts)', fontsize=20)
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.colorbar(label='dt $[\mu s]$')

    circ=plt.Circle((0,0), radius=CutLimits['RadialLimit'], color='mediumturquoise', fill=False)
    ax.add_patch(circ)

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'Radial_' + FilelistName + '.png')
    return fig


def PlotRadiusFlat(data, sPathToFigureSaveDirectory):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.hist(data.r**2, bins=np.linspace(0, 1000, 100))


    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title('FilelistName' + ': distribution in r$^2$ (all cuts)')
    plt.xlabel('$r^2$ [cm$^2$]')
    plt.ylabel('Counts')
    
    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'Radial_Hist_' + FilelistName + '.png')
    return fig


def PlotDepth(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))

    data = hax.cuts.range_selection(data, 's1', (CutLimits['LowLimitS1'], CutLimits['UpLimitS1']))
    data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))

    plt.scatter(data.r**2, data.z,
                c=data.s1,
                marker='.', s=6, color='k', linewidth=0, cmap='viridis')

    plt.hlines([CutLimits['LowBoundZ'], CutLimits['UpBoundZ']], 0,
               CutLimits['RadialLimit']**2, linestyle='solid', color='mediumturquoise')
    plt.vlines([CutLimits['RadialLimit']**2], CutLimits['LowBoundZ'],
               CutLimits['UpBoundZ'], linestyle='solid', color='mediumturquoise')
    plt.xlim(0,2500)
    plt.ylim(-105,5)
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': R vs. Z (S1 & S2 cuts)', fontsize=20)
    plt.xlabel('$r^2$ $[cm^2]$')
    plt.ylabel('z [cm]')
    plt.colorbar(label='S1 [PE]')

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'Depth_' + FilelistName + '.png')
    return fig


def PlotDepthFlat(data, sPathToFigureSaveDirectory):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1,1,1)
    ax.hist(data.z, bins=np.linspace(-100, 0, 100))

    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title('FilelistName' + ': distribution in z (all cuts)')
    plt.xlabel('z [cm]')
    plt.ylabel('Counts')
    
    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'Depth_Hist_' + FilelistName + '.png')
    return fig



def PlotS1S1AsymWithCuts(data, CutLimits, sPathToFigureSaveDirectory=''):
    fig = plt.figure(figsize=(10,8))

    data = hax.cuts.range_selection(data, 's2', (CutLimits['LowLimitS2'], CutLimits['UpLimitS2']))
    data = hax.cuts.range_selection(data, 's2_aft', (CutLimits['LowLimitAsymS2'], CutLimits['UpLimitAsymS2']))
    data = hax.cuts.selection(data, data.r < CutLimits['RadialLimit'], 'r2 < %i' %(CutLimits['RadialLimit'])**2)

#    Xs2asym = EvalCut(data, cuts, 'Xs2asym')
#    Xr = EvalCut(data, cuts, 'Xr')
#    Xz = EvalCut(data, cuts, 'Xz')
#    Xcuts = Xs2 & Xs2asym & Xr & Xz

    plt.scatter(data.s1, data.s1_aft,
                c=data.z, marker='.', edgecolors=None,
                vmin=-100, vmax=0, linewidths=0, s=6, cmap='viridis')

    plt.xlim(0,100000)
    plt.ylim(0,0.4)
    FilelistName = GetFilelistNameFromDirectory(sPathToFigureSaveDirectory)
    plt.title(FilelistName + ': S1 vs. S1 Asymmetry (Xs2, Xs2asym, Xr)', fontsize=20)
    plt.xlabel('S1 [PE]')
    plt.ylabel('S1 Area Fraction Top')
    plt.colorbar(label='z [cm]')

    plt.tight_layout()
    if len(sPathToFigureSaveDirectory) > 1:
        plt.savefig(sPathToFigureSaveDirectory + 'S1S1AsymWithCuts_' + FilelistName + '.png')
    return fig
