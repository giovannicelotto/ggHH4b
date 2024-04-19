import numpy as np
import uproot
import glob
import pandas as pd
import sys
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")

def main():
    print("Inizia la main")
    pathSignal = "/t3home/gcelotto/ggHH4b/flatData/myfileSig.parquet"
    pathBackground = "/t3home/gcelotto/ggHH4b/flatData/myfileBkg.parquet"
    dfSignal = pd.read_parquet(pathSignal)
    dfBackground = pd.read_parquet(pathBackground)
    print(dfSignal, dfBackground)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    
    ax[0].hist(dfSignal.mHH, histtype=u'step', color='blue', bins=np.linspace(0, 1500, 50))
    ax[0].hist(dfBackground.mHH, histtype=u'step', color='red', bins=np.linspace(0, 1500, 50))

    ax[1].hist(dfSignal.mH1, histtype=u'step', color='blue', bins=np.linspace(0, 200, 50))
    ax[1].hist(dfBackground.mH1, histtype=u'step', color='red', bins=np.linspace(0, 200, 50))

    ax[2].hist(dfSignal.mH2, histtype=u'step', color='blue', bins=np.linspace(0, 200, 50))
    ax[2].hist(dfBackground.mH2, histtype=u'step', color='red', bins=np.linspace(0, 200, 50))
    outName = "/t3home/gcelotto/ggHH4b/plots/3features.png"
    fig.savefig(outName, bbox_inches='tight')



    return 0

if __name__=="__main__":
    main()