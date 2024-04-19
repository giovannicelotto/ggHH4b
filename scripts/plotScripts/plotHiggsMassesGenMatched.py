import numpy as np
import uproot
import glob
import sys
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import linear_sum_assignment
from mpl_toolkits.axes_grid1 import make_axes_locatable
hep.style.use("CMS")
def main(ptThreshold):
    path = "/t3home/gcelotto/ggHH4b/flatData/df_mH1mH2_genMatched.parquet"
    df = pd.read_parquet(path)
    valRegion = df.lowestJetPt > ptThreshold
    print(np.sum(valRegion)/len(df))

    print(df.mH1.mean(), df.mH2.mean())
    print(df[valRegion].mH1.mean(), df[valRegion].mH2.mean())


    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    bins= np.linspace(0, 200, 30)
    ch1 = np.histogram(df.mH1, bins=bins)[0]
    ch2 = np.histogram(df.mH2, bins=bins)[0]
    ax[0].hist(bins[:-1], bins=bins, weights=ch1, label='All Events', alpha=0.6)
    ax[1].hist(bins[:-1], bins=bins, weights=ch2, label='All Events', alpha=0.6)
    ax[0].set_xlabel("m(H1) [GeV]")
    ax[1].set_xlabel("m(H2) [GeV]")

    ax[0].hist(df.mH1[valRegion], bins=bins, alpha=0.7, label='Events Jet Pt > %d GeV'%ptThreshold)
    ax[1].hist(df.mH2[valRegion], bins=bins, alpha=0.7, label='Events Jet Pt > %d GeV'%ptThreshold)

    for axs in ax:
        axs.set_ylim(axs.get_ylim())
        axs.vlines(x=125, ymin=0, ymax=axs.get_ylim()[1], color='black', linestyle='dotted', label='125 GeV')
        axs.vlines(x=120, ymin=0, ymax=axs.get_ylim()[1], color='red', linestyle='dotted', label='120 GeV')
        axs.legend(loc='upper left', fontsize=16)
        axs.set_ylabel("Events")
    ax[0].text(x=0, y=0.5, s="Mean : %.1f (%.1f)"%(df.mH1.mean(), df[valRegion].mH1.mean()), transform=ax[0].transAxes, fontsize=14)
    ax[1].text(x=0, y=0.5, s="Mean : %.1f (%.1f)"%(df.mH2.mean(), df[valRegion].mH2.mean()), transform=ax[1].transAxes, fontsize=14)
    
    outPath="/t3home/gcelotto/ggHH4b/plots/HHPeakGenMatched_pt%d.png"%ptThreshold
    print("saving in ", outPath)
    fig.savefig(outPath, bbox_inches='tight')

    return


if __name__ == "__main__":
    ptThreshold = int(sys.argv[1]) if len(sys.argv)>1 else None
    if ptThreshold == None:
        sys.exit("Specify ptThreshold")
    main(ptThreshold=ptThreshold)