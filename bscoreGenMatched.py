import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:10]
    print("%d files to be used" %len(fileNames))
    matchedEvents, totalEntries = 0, 0
    list_btagMatched=[]
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
            
            m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            if np.sum(m)!=4:
                continue
            matchedEvents=matchedEvents+1
            
            h1_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]==GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
            h2_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]!=GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
# np.arange(nJet)[m] are the recojets matched
            hh_RecoJetsIdx = np.arange(nJet)[m]
            h1_RecoJetsIdx = [idx for idx in hh_RecoJetsIdx if GenJet_partonMotherIdx[Jet_genJetIdx[idx]]==GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
            h2_RecoJetsIdx = [idx for idx in hh_RecoJetsIdx if GenJet_partonMotherIdx[Jet_genJetIdx[idx]]!=GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
            list_btagMatched.append(Jet_btagDeepFlavB[hh_RecoJetsIdx])
            

            
    
    fig, ax = plt.subplots(1, 1)
    bins=np.linspace(0, 1, 20)
    c = np.histogram(list_btagMatched, bins=bins)[0]
    c=c/np.sum(c)
    ax.hist(bins[:-1], bins=bins, weights=c)
    outName = "/t3home/gcelotto/ggHH4b/hh_btagMatched.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Savin in", outName)
    print(matchedEvents, matchedEvents/totalEntries)
    return

if __name__ == "__main__":
    main()