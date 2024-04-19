'''
script to plot the effect of requiring a triggering muon in the 4 jets from higgs
'''
import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
hep.style.use("CMS")


def trigMuon():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:20]
    print("%d files to be used" %len(fileNames))
    ntrigMu = []
    trigJet = []
    nTotalJet = []
    for fileName in fileNames:
        print("\nOpening %d/%d"%(fileNames.index(fileName)+1, len(fileNames)), fileName)
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries   
        for ev in  range(maxEntries):
            if (ev%(int(maxEntries/20))==0):
                sys.stdout.write('\r')
                sys.stdout.write("%d%%"%(ev/maxEntries*100))
                sys.stdout.flush()
            Jet_muonIdx1        = branches["Jet_muonIdx1"][ev]
            Jet_muonIdx2        = branches["Jet_muonIdx2"][ev]
            Muon_isTriggering   = branches["Muon_isTriggering"][ev]
            Jet_eta             = branches["Jet_eta"][ev]
            nJet                = branches["nJet"][ev]
            nTrigJetEvent       = 0
            nTotalJet.append(nJet)

            for jetIdx in range(nJet):
                if Jet_muonIdx1[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx1[jetIdx]]):
                        nTrigJetEvent=nTrigJetEvent+1
                        continue
                if Jet_muonIdx2[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx2[jetIdx]]):
                        nTrigJetEvent=nTrigJetEvent+1
                        continue
            trigJet.append(nTrigJetEvent)
            ntrigMu.append(np.sum(Muon_isTriggering))

    
    # Create a 2D histogram
    hist, xedges, yedges = np.histogram2d(ntrigMu, trigJet,  bins=(np.arange(5)-0.5, np.arange(5)-0.5))

    fig, ax = plt.subplots(1, 1)
    total_count = np.sum(hist)
    print(hist/total_count*100)
    im = ax.imshow(hist/total_count, cmap='RdYlBu', origin='lower')# extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]
    for i in range(4):
        for j in range(4):
            text = f'{hist[j, i] / total_count:.1%}'  # Calculate percentage
            if (j ==0) | (i>j):
                ax.text(xedges[i] + 0.5, yedges[j] + 0.5, "--", color='black', ha='center', va='center')
            else:
                ax.text(xedges[i] + 0.5, yedges[j] + 0.5, text, color='black', ha='center', va='center')

    cbar = fig.colorbar(im, ax=ax, label='Percentage', shrink=0.8)  # Add color bar
    ax.set_xlabel('nTrigJets')
    ax.set_ylabel('nTrigMuon')
    ax.set_xticks(range(4))
    ax.set_xticklabels([str(i) for i in range(4)])
    ax.set_yticks(range(4))
    ax.set_yticklabels([str(i) for i in range(4)])

    # Manually set the color of a specific cell (e.g., cell at row 2, column 3) to red
    for col in range(4):
        for row in range(col):
            rect = plt.Rectangle((xedges[col], yedges[row]), xedges[col + 1] - xedges[col], yedges[row + 1] - yedges[row], fill=True, color='white', alpha=0.8)
            ax.add_patch(rect)
    rect = plt.Rectangle((xedges[0], yedges[0]), xedges[col + 1] - xedges[col], yedges[row + 1] - yedges[row], fill=True, color='white', alpha=0.8)
    ax.add_patch(rect)

    fig.savefig("/t3home/gcelotto/ggHH4b/matrixTrigMuonsvsTrigJets.png")



# Jets multiplicity
    fig, ax = plt.subplots(1, 1)
    bins = np.arange(38)
    c = np.histogram(nTotalJet, bins=bins)[0]
    c=c/np.sum(c)*100
    ax.hist(bins[:-1], bins=bins, weights=c, color='blue', histtype=u'step')
    ax.set_ylabel("Percentage [%]")
    ax.set_xlabel("N$_\mathrm{jets}$")
    ax.text(x=0.95, y=0.9, s="Mean    %.2f"%np.mean(nTotalJet), ha='right', transform=ax.transAxes)
    ax.text(x=0.95, y=0.84, s="Std Dev    %.2f"%np.std(nTotalJet), ha='right', transform=ax.transAxes)
    outName = "/t3home/gcelotto/ggHH4b/nJets.png"
    fig.savefig(outName, bbox_inches='tight')

# Triggering Muon multiplicity
    fig, ax = plt.subplots(1, 1)
    bins = np.arange(6)
    c = np.histogram(ntrigMu, bins=bins)[0]
    c=c/np.sum(c)*100
    ax.hist(bins[:-1], bins=bins-0.5, weights=c, color='blue', histtype=u'step')
    ax.set_ylabel("Percentage [%]")
    ax.set_xlabel("N$_\mathrm{trig \mu}$")
    ax.text(x=0.95, y=0.9, s="Mean    %.2f"%np.mean(ntrigMu), ha='right', transform=ax.transAxes)
    ax.text(x=0.95, y=0.84, s="Std Dev    %.2f"%np.std(ntrigMu), ha='right', transform=ax.transAxes)
    outName = "/t3home/gcelotto/ggHH4b/nTrigMuon.png"
    fig.savefig(outName, bbox_inches='tight')

# Triggering Jets multiplicity
    fig, ax = plt.subplots(1, 1)
    bins = np.arange(6)
    c = np.histogram(trigJet, bins=bins)[0]
    c=c/np.sum(c)*100
    ax.hist(bins[:-1], bins=bins-0.5, weights=c, color='blue', histtype=u'step')
    ax.set_ylabel("Percentage [%]")
    ax.set_xlabel("N$_\mathrm{trig jets}$")
    ax.text(x=0.95, y=0.9, s="Mean    %.2f"%np.mean(trigJet), ha='right', transform=ax.transAxes)
    ax.text(x=0.95, y=0.84, s="Std Dev    %.2f"%np.std(trigJet), ha='right', transform=ax.transAxes)
    outName = "/t3home/gcelotto/ggHH4b/nTrigJets.png"
    fig.savefig(outName, bbox_inches='tight')
    return 0


if __name__ == "__main__":
    trigMuon()