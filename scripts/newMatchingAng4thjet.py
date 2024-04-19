import numpy as np
import uproot
import glob
import sys
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import linear_sum_assignment
from mpl_toolkits.axes_grid1 import make_axes_locatable
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:10]
    print("%d files to be used" %len(fileNames))
    matchedEvents, totalEntries = 0, 0, 
    completeMatching = 0
    partialMatching = 0
    partialNotWorking = 0
    angDiff = []
    firstJetsDiff = []
    for fileName in fileNames:
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries 
        totalEntries = totalEntries + maxEntries
        print("Entries : %d" %maxEntries)

        for ev in  range(maxEntries):
            
            GenJet_partonFlavour        = branches["GenJet_partonFlavour"][ev]
            GenJet_partonMotherIdx      = branches["GenJet_partonMotherIdx"][ev]
            GenJet_partonMotherPdgId    = branches["GenJet_partonMotherPdgId"][ev]

            GenPart_statusFlags         = branches["GenPart_statusFlags"][ev]
            GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
            GenPart_pdgId               = branches["GenPart_pdgId"][ev]
            GenPart_eta                 = branches["GenPart_eta"][ev]
            GenPart_phi                 = branches["GenPart_phi"][ev]
            GenPart_mass                = branches["GenPart_mass"][ev]
        # Reco Jets
            nJet                        = branches["nJet"][ev]
            Jet_eta                     = branches["Jet_eta"][ev]
            Jet_pt                      = branches["Jet_pt"][ev]
            Jet_phi                     = branches["Jet_phi"][ev]
            Jet_mass                    = branches["Jet_mass"][ev]
            Jet_bReg2018                 = branches["Jet_bReg2018"][ev]
            Jet_genJetIdx               = branches["Jet_genJetIdx"][ev]
            Jet_btagDeepFlavB           = branches["Jet_btagDeepFlavB"][ev]
            GenJet_pt                   = branches["GenJet_pt"][ev]
            GenJet_eta                  = branches["GenJet_eta"][ev]
            GenJet_phi                  = branches["GenJet_phi"][ev]
            GenJet_mass                 = branches["GenJet_mass"][ev]
            Jet_muonIdx1                = branches["Jet_muonIdx1"][ev]
            Jet_muonIdx2                = branches["Jet_muonIdx2"][ev]
            Muon_isTriggering           = branches["Muon_isTriggering"][ev]
            Jet_qgl                     = branches["Jet_qgl"][ev]
            
            # limit the data to events where 4 jets are gen matched
            m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            if np.sum(m)==4:
                completeMatching = completeMatching+1
                pass
            elif np.sum(m)==3:
                #print(list(zip(Jet_eta[m], Jet_phi[m])))
                
                newM = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5)
                
                if np.sum(newM)==4:
                    # check if the fourth jet found is really the one missing from the quark
                    mQuark = (GenPart_pdgId[GenPart_genPartIdxMother]==25) & ((GenPart_statusFlags[GenPart_genPartIdxMother] & 8192)!=0)
                    Jets = np.array(list(zip(Jet_eta[newM], Jet_phi[newM])))
                    Quarks = np.array(list(zip(GenPart_eta[mQuark], GenPart_phi[mQuark])))

                # First control
                    newJetMask = newM ^ m # or esclusivo tra i due, true nel quarto jet aggiunto
                    newJetAdded = np.arange(4)[Jet_eta[newM] == Jet_eta[newJetMask]][0]  # index tra 0 e 3 del quarto jet aggiunto
                    
                    distances = np.linalg.norm(Quarks[:, np.newaxis, :] - Jets, axis=2)
                    #print("DISTANCES:\n")
                    #print(distances)

                    # Use the Hungarian algorithm to find the optimal assignment
                    QuarkIdx, JetIdx = linear_sum_assignment(distances)

                    # Create a dictionary to store the matching pairs
                    matching_pairs = {i: j for i, j in zip(JetIdx, QuarkIdx)}
                    #print("New Jet added : ", newJetAdded) 
                    #print("New Quark added : ", matching_pairs[newJetAdded]) 
                    #print(Jets)
                    #print(Quarks)
                    #print(matching_pairs)
                    #print("Jet eta : ",   Jet_eta[newM][newJetAdded])
                    #print("Jet phi : ", Jet_phi[newM][newJetAdded])
                    #print("Quark eta : ",   GenPart_eta[mQuark][matching_pairs[newJetAdded]])
                    #print("Quark phi : ", GenPart_phi[mQuark][matching_pairs[newJetAdded]])
                    #input("Next")
                    angDiff.append([Jet_eta[newM][newJetAdded] - GenPart_eta[mQuark][matching_pairs[newJetAdded]], Jet_phi[newM][newJetAdded] - GenPart_phi[mQuark][matching_pairs[newJetAdded]]])
                    for i in np.arange(4):
                        if i == newJetAdded:
                            continue
                        else:
                            firstJetsDiff.append([Jet_eta[newM][i] - GenPart_eta[mQuark][matching_pairs[i]], Jet_phi[newM][i] - GenPart_phi[mQuark][matching_pairs[i]]])



                    partialMatching = partialMatching +1
                    m = newM
                elif np.sum(newM)!=4:
                    partialNotWorking = partialNotWorking + 1

            


    

    print("Partial Matching : ", partialMatching/totalEntries)
    print("Complete Matching : ", completeMatching/totalEntries)
    print("Partial not working : ", partialNotWorking/totalEntries)


    etaDiff, phiDiff = np.array(angDiff)[:,0], np.array(angDiff)[:,1]
    x_bins, y_bins = np.linspace(-2, 2, 30), np.linspace(-np.pi/2, np.pi/2, 30)
    fig, ax_main = plt.subplots(figsize=(8, 8))
    divider = make_axes_locatable(ax_main)
    ax_top = divider.append_axes("top", 1.2, pad=0.2, sharex=ax_main)
    ax_right = divider.append_axes("right", 1.2, pad=0.2, sharey=ax_main)
    ax_top.set_title("4th matched jet")

    # Plot the 2D histogram in the main axes
    hist, x_edges, y_edges = np.histogram2d(x=etaDiff, y=phiDiff, bins=[x_bins, y_bins])
    ax_main.imshow(hist.T, origin='lower', extent=(x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()), aspect='auto', cmap='Blues')
    ax_main.set_xlabel("Jet eta - Quark eta")
    ax_main.set_ylabel("Jet phi - Quark phi")

    # Plot the marginalized histogram on top
    ax_top.hist(etaDiff, bins=x_bins, color='lightblue', edgecolor='black')
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_yticks([])
    ax_top.xaxis.tick_top()

    # Plot the marginalized histogram on the right
    ax_right.hist(phiDiff, bins=y_bins, color='lightblue', edgecolor='black', orientation='horizontal')#lightcoral
    ax_right.set_ylim(ax_main.get_ylim())
    ax_right.set_xticks([])
    ax_right.yaxis.tick_right()
    outName = "/t3home/gcelotto/ggHH4b/plots/4thJetVsQuark.png"
    fig.savefig(outName, bbox_inches='tight')
    print("saving in ", outName)




    # same plot for the other 3 jets
    etaDiff, phiDiff = np.array(firstJetsDiff)[:,0], np.array(firstJetsDiff)[:,1]
    x_bins, y_bins = np.linspace(-2, 2, 30), np.linspace(-np.pi/2, np.pi/2, 30)
    fig, ax_main = plt.subplots(figsize=(8, 8))
    divider = make_axes_locatable(ax_main)
    ax_top = divider.append_axes("top", 1.2, pad=0.2, sharex=ax_main)
    ax_right = divider.append_axes("right", 1.2, pad=0.2, sharey=ax_main)
    ax_top.set_title("First 3 matched jets ")

    # Plot the 2D histogram in the main axes
    hist, x_edges, y_edges = np.histogram2d(x=etaDiff, y=phiDiff, bins=[x_bins, y_bins])
    ax_main.imshow(hist.T, origin='lower', extent=(x_bins.min(), x_bins.max(), y_bins.min(), y_bins.max()), aspect='auto', cmap='Blues')
    ax_main.set_xlabel("Jet eta - Quark eta")
    ax_main.set_ylabel("Jet phi - Quark phi")

    # Plot the marginalized histogram on top
    ax_top.hist(etaDiff, bins=x_bins, color='lightblue', edgecolor='black')
    ax_top.set_xlim(ax_main.get_xlim())
    ax_top.set_yticks([])
    ax_top.xaxis.tick_top()

    # Plot the marginalized histogram on the right
    ax_right.hist(phiDiff, bins=y_bins, color='lightblue', edgecolor='black', orientation='horizontal')#lightcoral
    ax_right.set_ylim(ax_main.get_ylim())
    ax_right.set_xticks([])
    ax_right.yaxis.tick_right()
    outName = "/t3home/gcelotto/ggHH4b/plots/3jetsVsQuark.png"
    fig.savefig(outName, bbox_inches='tight')
    print("saving in ", outName)

    print("len angdiff", len(angDiff))
    print("len 3jets ang diff", len(firstJetsDiff))
    return

if __name__ == "__main__":
    main()