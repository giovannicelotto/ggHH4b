import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:20]
    print("%d files to be used" %len(fileNames))

    matchedEvents, totalEntries = 0, 0
    # numerator and denominator for the ratio
    trigJetFromHiggs = 0
    trigJet = 0

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
            Jet_genJetIdx               = branches["Jet_genJetIdx"][ev]
            nJet                        = branches["nJet"][ev]
            Jet_muonIdx1                = branches["Jet_muonIdx1"][ev]
            Jet_muonIdx2                = branches["Jet_muonIdx2"][ev]
            Muon_isTriggering           = branches["Muon_isTriggering"][ev]

            # mask for recoJets to matched recoJets
            m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            if np.sum(m)!=4:
                continue

            matchedEvents=matchedEvents+1

            for jetIdx in range(nJet):
                # loop over all the jes.
                # find trig jets.

                # check if trig jets in matched jets
                trigJetFilled = False
                if Jet_muonIdx1[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx1[jetIdx]]):
                        trigJet = trigJet + 1
                        trigJetFilled = True
                        if jetIdx in np.arange(nJet)[m]:
                            trigJetFromHiggs = trigJetFromHiggs + 1
                            continue
                if Jet_muonIdx2[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx2[jetIdx]]):
                        if not trigJetFilled:
                            trigJet = trigJet + 1
                        if jetIdx in np.arange(nJet)[m]:
                            trigJetFromHiggs = trigJetFromHiggs + 1
                            continue
    print(trigJetFromHiggs/trigJet)

if __name__ == "__main__":
    main()
            


