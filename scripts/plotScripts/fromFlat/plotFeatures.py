import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import mplhep as hep
hep.style.use("CMS")
'''

nFiles = int(sys.argv[1]) if len(sys.argv)>1 else None

'''
def plotFeatures(signal, realData, outName):
    
    

    nRow, nCol = 12, 8
    fig, ax = plt.subplots(nRow, nCol, figsize=(30, 40), constrained_layout=True)
    fig.align_ylabels(ax[:,0])
    xlims = np.array([
        [0, 400, 30], [-5, 5, 30], [0, np.pi, 30], [0, 1000, 30],
        [0, 400, 30], [-5, 5, 30], [-np.pi, np.pi, 30], [0, 200, 30],
        [0, 300, 30], [-5, 5, 30], [0, np.pi, 30], [0, 200, 30],
        [0, 250, 30], [-5, 5, 30], [0, np.pi, 30], [0, 30, 30],
        [0, 250, 30], [-5, 5, 30], [0, np.pi, 30], [0, 30, 30],
        [0, 250, 30], [-5, 5, 30], [0, np.pi, 30], [0, 30, 30],
        [0, 100, 30], [-5, 5, 30], [0, np.pi, 30], [0, 30, 30],

        [0, 1, 10],     [0, 1, 10], [0, 1, 10],         [0, 1, 10],
        [0, 1, 10],     [0, 1, 10], [0, 1, 10],         [0, 1, 10],
        [0, 50, 10],    [0, 50, 10],[0, 50, 10],        [0, 50, 10],
        [0, 4, 5],     [0, 4, 5], [0, 4, 5],         [0, 4, 5],
        [0, 4, 5],     [0, 4, 5], [0, 4, 5],         [0, 4, 5],

                        [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
                        [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
                        [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, np.pi, 20], [0, 5, 20], [0, 6, 20], [0, np.pi/2, 20],
        [0, 900, 30],
        [0, 100, 30], [-5, 5, 30], [-np.pi, np.pi, 30], [0.102, 0.11, 30],
    ])
    labels = signal.columns.tolist()
    for i in range(nRow):
        fig.align_xlabels(ax[i,:])
        for j in range(nCol):
            fig.align_ylabels(ax[:,j])
            if i*nCol+j>=len(signal.columns):
                break
            bins = np.linspace(xlims[i*nCol+j,0], xlims[i*nCol+j,1], int(xlims[i*nCol+j,2]) + 1)
            
            countsSignal = np.zeros(len(bins) - 1)
            countsSignal = np.histogram(np.clip(signal.iloc[:,i*nCol+j], bins[0], bins[-1]), bins=bins)[0]
            countsSignalErr = np.sqrt(countsSignal)
            
            countsBkg = np.zeros(len(bins)-1)    
            countsBkg = np.histogram(np.clip(realData.iloc[:,i*nCol+j], bins[0], bins[-1]), bins=bins)[0]
            countsBkgErr = np.sqrt(countsBkg)

            # Normalize the counts to 1 so also the errors undergo the same operation. Do first the errors, otherwise you lose the info on the signal
            countsSignalErr = countsSignalErr/np.sum(countsSignal)
            countsSignal = countsSignal/np.sum(countsSignal)
            countsBkgErr=countsBkgErr/np.sum(countsBkg)
            countsBkg=countsBkg/np.sum(countsBkg)

            ax[i, j].hist(bins[:-1], bins=bins, weights=countsSignal,   label='Signal', histtype=u'step',  color='blue', )
            ax[i, j].hist(bins[:-1], bins=bins, weights=countsBkg,      label='Background', histtype=u'step' , color='red',)

            ax[i, j].set_xlabel(labels[i*nCol+j], fontsize=18)
            ax[i, j].set_xlim(bins[0], bins[-1])
            ax[i, j].set_ylabel("Probability", fontsize=18)

            ax[i, j].legend(fontsize=18)
            ax[i, j].tick_params(which='major', length=8)
            ax[i, j].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i, j].tick_params(which='minor', length=4)
            
    
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in %s"%outName)
    plt.close('all')

    return

if __name__=="__main__":


    nFiles = int(sys.argv[1]) if len(sys.argv)>1 else None
    if nFiles == None:
        sys.exit("Specify nFiles (int)")
    realDataPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/flatDataForggHH4b"
    signalPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData"
    signalFileNames = glob.glob(signalPath+"/*.parquet")
    realDataFileNames = glob.glob(realDataPath+"/*.parquet")[:nFiles]

    signal = pd.read_parquet(signalFileNames)
    realData = pd.read_parquet(realDataFileNames)
    outName = "/t3home/gcelotto/ggHH4b/plots/features/features.png"
    plotFeatures(signal=signal, realData=realData, outName=outName)