import numpy as np
import glob
import pandas as pd
import csv
import sys
sys.path.append('/t3home/gcelotto/ggHH4b/scripts')
from helpers import loadDask, loadParquet, getXSectionBR
from plotFeatures import plotFeatures 
import dask.dataframe as dd 
import dask
import mplhep as hep
hep.style.use("CMS")
dask.config.set({'distributed.worker.memory.target': 0.8})
dask.config.set(scheduler='threads', num_workers=4)

def cutDaskDataFrame(signal, realData, featureName, min, max):
    maskSignal = (signal.h1_mass)>-9999 # always true
    if min is not None:
        maskSignal = (maskSignal) & (signal[featureName] > min)
    if max is not None:
        maskSignal = (maskSignal) & (signal[featureName] < max)

    maskData = (realData.h1_mass)>-9999 # always true
    if min is not None:
        maskData = (maskData) & (realData[featureName] > min)
    if max is not None:
        maskData = (maskData) & (realData[featureName] < max)

    return signal[maskSignal], realData[maskData]

def is_point_inside_ellipse(x, y, center_x, center_y, a, b):
        distance = ((x - center_x) / a)**2 + ((y - center_y) / b)**2
        return distance <= 1

def cutAndCount(nRealDataFiles):
    realDataPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/flatDataForggHH4b"
    signalPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData"
    nRealDataFiles = len(glob.glob(realDataPath+"/*.parquet")) if nRealDataFiles==-1 else nRealDataFiles
    currentLumi = 0.774*nRealDataFiles/1017
    N_SignalMini = np.load("/t3home/gcelotto/ggHH4b/outputs/N_mini.npy")
    correctionSignal = 1/N_SignalMini*getXSectionBR()*currentLumi
    signal, realData = loadDask(signalPath=signalPath, realDataPath=realDataPath, nSignalFiles=-1, nRealDataFiles=nRealDataFiles)
    startingSignal, startingBkg = len(signal), len(realData)
    print("initial signal : ", startingSignal)
    print("initial bkg : ", startingBkg)
    
    # apply cuts
    with open('/t3home/gcelotto/ggHH4b/outputs/ellipseCoord.csv', 'r', newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        for row in csv_reader:
            
            if row and row[0].startswith('#'):
                pass
            else:
                x_max, y_max, sigma_x, sigma_y = [float(num) for num in row]
    signalMask = is_point_inside_ellipse(signal.h1_mass, signal.h2_mass, x_max, y_max, sigma_x, sigma_y)
    realDataMask = is_point_inside_ellipse(realData.h1_mass, realData.h2_mass, x_max, y_max, sigma_x, sigma_y)
    signal = signal[signalMask]
    realData = realData[realDataMask]


    signal, realData = cutDaskDataFrame(signal,realData, 'j1_eta', -2.5, 2.5)
    signal, realData = cutDaskDataFrame(signal,realData, 'j2_eta', -2.5, 2.5)
    signal, realData = cutDaskDataFrame(signal,realData, 'j3_eta', -2.5, 2.5)
    signal, realData = cutDaskDataFrame(signal,realData, 'j4_eta', -2.5, 2.5)
    signal, realData = cutDaskDataFrame(signal,realData, 'j1_pt', 20, None)
    signal, realData = cutDaskDataFrame(signal,realData, 'j2_pt', 20, None)
    signal, realData = cutDaskDataFrame(signal,realData, 'j3_pt', 20, None)
    signal, realData = cutDaskDataFrame(signal,realData, 'j4_pt', 20, None)
    #signal, realData = cutDaskDataFrame(signal,realData, 'h1_pt', 100, None)
    #signal, realData = cutDaskDataFrame(signal,realData, 'h2_pt', 70, None)
    #signal, realData = cutDaskDataFrame(signal,realData, 'h1_pt', 300, None)
    realData = realData.compute()
    signal = signal.compute()

    

    
    finalSignal, finalBkg = len(signal), len(realData)
    print("final signal : ",finalSignal, finalSignal*correctionSignal,  finalSignal/startingSignal*100, " %")
    print("final signal full lumi: ",finalSignal*correctionSignal*41.6/currentLumi)
    print("final bkg : ",finalBkg,  finalBkg/startingBkg*100, " %")
    print("S/B : ", finalSignal*correctionSignal/finalBkg)
    print("S/sqrt(B) : ", finalSignal*correctionSignal/np.sqrt(finalBkg))
    print("S/sqrt(B) full lumi: ", finalSignal*correctionSignal/np.sqrt(finalBkg)*np.sqrt(41.6/currentLumi))

    outName = "/t3home/gcelotto/ggHH4b/plots/features/featuresCut.png"
    plotFeatures(signal=signal, realData=realData, outName=outName)
    
  


    return



if __name__ =="__main__":
    nRealDataFiles = int(sys.argv[1]) if len(sys.argv)>1 else None
    if nRealDataFiles == None:
        sys.exit("Specify nRealDataFiles (int)")
    cutAndCount(nRealDataFiles=nRealDataFiles)