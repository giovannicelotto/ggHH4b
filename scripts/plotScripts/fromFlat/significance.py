import matplotlib.pyplot as plt
import numpy as np
import glob
import csv
from matplotlib.ticker import AutoMinorLocator, LogLocator
import matplotlib.patches as patches
from scipy.optimize import curve_fit
from scipy.stats import norm, crystalball
from scipy.integrate import quad
import sys
sys.path.append('/t3home/gcelotto/ggHH4b/scripts')
from helpers import loadDask, loadParquet, getXSectionBR
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm
import mplhep as hep
hep.style.use("CMS")
def is_point_inside_ellipse(x, y, center_x, center_y, a, b):
    distance = ((x - center_x) / a)**2 + ((y - center_y) / b)**2
    return distance <= 1

def add_text_to_hist(ax, bins, h):
    for i in range(len(bins) - 1):
        for j in range(len(bins) - 1):
            count = h[0][i, j]
            if count > 0:
                ax.text(
                    bins[i] + (bins[i + 1] - bins[i]) / 2,    bins[j] + (bins[j + 1] - bins[j]) / 2,
                    f'{count:.0f}', color='black',
                    ha='center', va='center',
                    fontsize=8
                )

def plotmH1(realFiles=1):
    signalPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData"
    realDataPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/flatDataForggHH4b"
    signal, background = loadParquet(signalPath, realDataPath, nSignalFiles=-1, nRealDataFiles=realFiles, columns=['h1_mass', 'h2_mass'])
    print(len(background), " events for Parking")
    print(len(signal), " events for ggHH")
    realFiles = len(glob.glob(realDataPath+"/*.parquet")) if realFiles==-1 else realFiles
    currentLumi = 0.774*realFiles/1017
    N_SignalMini = np.load("/t3home/gcelotto/ggHH4b/outputs/N_mini.npy")
    correctionSignal = 1/N_SignalMini*getXSectionBR()*currentLumi
    print("Current Lumi ", currentLumi)
    print("Correction signal", correctionSignal)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    bins = np.linspace(0, 220, 40)
# mH1
    visibilityCorrection = 10**7
    cSig = np.histogram(np.clip(signal.h1_mass, bins[0], bins[-1]), bins=bins)[0]*correctionSignal*visibilityCorrection
    cBkg = np.histogram(np.clip(background.h1_mass, bins[0], bins[-1]), bins=bins)[0]
    ax[0].hist(bins[:-1], bins=bins, weights=cSig, color='blue', histtype=u'step', label=r'ggHH4b $\times 10^{7}$')
    ax[0].hist(bins[:-1], bins=bins, weights=cBkg, color='red', histtype=u'step', label='background')
    ax[0].set_ylim(ax[0].get_ylim()[0], ax[0].get_ylim()[1])
    ax[0].legend(loc='upper right', fontsize=16)
# mH2
    cSig = np.histogram(np.clip(signal.h2_mass, bins[0], bins[-1]), bins=bins)[0]*correctionSignal*visibilityCorrection
    cBkg = np.histogram(np.clip(background.h2_mass, bins[0], bins[-1]), bins=bins)[0]
    ax[1].hist(bins[:-1], bins=bins, weights=cSig, color='blue', histtype=u'step', label=r'ggHH4b $\times 10^{7}$')
    ax[1].hist(bins[:-1], bins=bins, weights=cBkg, color='red', histtype=u'step', label='background')
    ax[1].set_ylim(ax[1].get_ylim()[0], ax[1].get_ylim()[1])
    ax[1].legend(loc='best', fontsize=16)
    outName = "/t3home/gcelotto/ggHH4b/plots/significance/mH1_significance.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)
    
    # heatmap for the signal
    fig, ax = plt.subplots(1, 1)
    bins=np.linspace(50, 300, 40)
    h = ax.hist2d(signal.h1_mass, signal.h2_mass, weights=np.ones(len(signal))*correctionSignal, bins=(bins, bins), cmap='viridis', norm=LogNorm())

    max_index = np.unravel_index(np.argmax(h[0]), h[0].shape)
    x_max = h[1][max_index[0]] + (bins[1]-bins[0])/2
    y_max = h[2][max_index[1]] + (bins[1]-bins[0])/2
    sigma_x = 23.4
    sigma_y = 30.86
    
    with open('/t3home/gcelotto/ggHH4b/outputs/ellipseCoord.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["#x_max, y_max, sigma_x, sigma_y"])
        csv_writer.writerow([x_max, y_max, sigma_x, sigma_y])
    np.save("/t3home/gcelotto/ggHH4b/outputs/xmax_ymax_sigmax_sigmay.npy", [x_max, y_max, sigma_x, sigma_y])
    print("Maximum in the 2D plane : ", x_max, y_max)
    
    ax.vlines(x_max,ymin=0, ymax=bins[-1], color='red')
    ax.hlines(y_max,xmin=0, xmax=bins[-1], color='red')
    ax.add_patch(Ellipse((x_max, y_max), 2 * sigma_x, 2 * sigma_y, color='red', alpha=0.3, label='Ellipse'))
    ax.set_xlabel('mH1 [GeV]')
    ax.set_ylabel('mH2 [GeV]')
    cbar = plt.colorbar(h[3], ax=ax)
    cbar.set_label('Counts')
    outName = "/t3home/gcelotto/ggHH4b/plots/significance/mH1mH2_heatmap_signal.png"
    #add_text_to_hist(ax, bins, h)
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)

    #heatmap for background
    fig, ax = plt.subplots(1, 1)
    h = ax.hist2d(background.h1_mass, background.h2_mass, bins=(bins, bins), cmap='viridis', norm=LogNorm())
    ax.vlines(x_max,ymin=0, ymax=bins[-1], color='red')
    ax.hlines(y_max,xmin=0, xmax=bins[-1], color='red')
    ax.add_patch(Ellipse((x_max, y_max), 2 * sigma_x, 2 * sigma_y, color='red', alpha=0.3, label='Ellipse'))
    ax.set_xlabel('mH1 [GeV]')
    ax.set_ylabel('mH2 [GeV]')
    cbar = plt.colorbar(h[3], ax=ax)
    cbar.set_label('Counts')
    outName = "/t3home/gcelotto/ggHH4b/plots/significance/mH1mH2_heatmap.png"
    #add_text_to_hist(ax, bins, h)

    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)

    

    print("Significance at lumi : ", currentLumi, " fb-1")
    #signal_SR = np.sum((signal.h1_mass>mean_fit_h1-2*sigma_fit_h1 ) & (signal.h1_mass<mean_fit_h1+2*sigma_fit_h1 ) & (signal.h2_mass>mean_fit_h2-2*sigma_fit_h2 ) & (signal.h2_mass<mean_fit_h2+2*sigma_fit_h2 ))*correctionSignal
    #data_SR = np.sum((background.h1_mass>mean_fit_h1-2*sigma_fit_h1 ) & (background.h1_mass<mean_fit_h1+2*sigma_fit_h1 ) & (background.h2_mass>mean_fit_h2-2*sigma_fit_h2 ) & (background.h2_mass<mean_fit_h2+2*sigma_fit_h2 ))
    signal_SR = np.sum(is_point_inside_ellipse(signal.h1_mass, signal.h2_mass, x_max, y_max, sigma_x, sigma_y))*correctionSignal
    data_SR = np.sum(is_point_inside_ellipse(background.h1_mass, background.h2_mass, x_max, y_max, sigma_x, sigma_y))
    print("signal_SR : ",np.sum(is_point_inside_ellipse(signal.h1_mass, signal.h2_mass, x_max, y_max, sigma_x, sigma_y)), signal_SR)
    print("data_SR : ",data_SR)
    print("signal_SR full lumi : ",signal_SR*41.6/currentLumi)
    print("data_SR full lumi: ",data_SR*41.6/currentLumi)
    print("S/B : ", signal_SR/data_SR)
    print("S/sqrtB : ", signal_SR/np.sqrt(data_SR))
    print("S/sqrt(B) full lumi : ",signal_SR/np.sqrt(data_SR)*np.sqrt(41.6/currentLumi))
    #print("Events in window mH1", np.sum((background.h1_mass>x_max-2*sigma_x ) & (background.h1_mass<x_max+2*sigma_x )))
    #print("Events in window mH2", np.sum((background.h2_mass>y_max-2*sigma_y ) & (background.h2_mass<y_max+2*sigma_y )))
    print("Total events in dataframe : ", len(background))

    
    #s_over_b_max, a_max, b_max=0, 0, 0
    #for a in np.linspace(20, 30, 30):
    #    for b in np.linspace(25, 35, 30):
    #        signal_SR = np.sum(is_point_inside_ellipse(signal.h1_mass, signal.h2_mass, x_max, y_max, a, b))
    #        data_SR = np.sum(is_point_inside_ellipse(background.h1_mass, background.h2_mass, x_max, y_max, a, b))
    #        s_over_b_current = signal_SR/np.sqrt(data_SR)
    #        if s_over_b_current>s_over_b_max:
    #            s_over_b_max = s_over_b_current
    #            a_max = a
    #            b_max= b
            
    #        print(a, b, signal_SR/data_SR)
    #print("MAX FOUND IN", a_max, b_max, s_over_b_max)
    

    

    return

if __name__ == "__main__":
    realFiles= int(sys.argv[1]) if len(sys.argv)>1 else -1
    plotmH1(realFiles=realFiles)