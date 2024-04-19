import numpy as np
import uproot
import glob
import sys
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:6]
    print("%d files to be used" %len(fileNames))
    correct, matchedEvents, totalEntries = 0, 0, 0
    
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
            Jet_muonIdx1                = branches["Jet_muonIdx1"][ev]
            Jet_muonIdx2                = branches["Jet_muonIdx2"][ev]
            Muon_isTriggering           = branches["Muon_isTriggering"][ev]
            Jet_qgl                     = branches["Jet_qgl"][ev]
            
            # limit the data to events where 4 jets are gen matched
            m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            if np.sum(m)==4:
                matchedEvents=matchedEvents+1
            
            elif np.sum(m)==3:

                newM = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5)
                
                if np.sum(newM)==4:
                    m=newM
                    matchedEvents=matchedEvents+1
                else:
                    continue
            else:
                continue
            
            

            # list of bools [nJet]
            IsTrigJet = []

            # list of scores based on pT (max), mass (max), |eta| (min), btag (max)
            scores = []
            jetsToCheck = np.min((12, nJet))
            for jetIdx in range(jetsToCheck):
                #jetTemp = ROOT.TLorentzVector(0.,0.,0.,0.)
                #jetTemp.SetPtEtaPhiM(Jet_pt[jetIdx], Jet_eta[jetIdx], Jet_phi[jetIdx], Jet_mass[jetIdx])
                #myRootJets.append(jetTemp)
                
                if Jet_muonIdx1[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx1[jetIdx]]):
                        IsTrigJet.append(100)
                        continue
                if Jet_muonIdx2[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx2[jetIdx]]):
                        IsTrigJet.append(100)
                    else:
                        IsTrigJet.append(0)
                else:
                    IsTrigJet.append(0)
            
            scores = (Jet_pt[:jetsToCheck] - np.min(Jet_pt[:jetsToCheck])) / (np.max(Jet_pt[:jetsToCheck]) - np.min(Jet_pt[:jetsToCheck]))
            scores = scores + (1- (abs(Jet_eta[:jetsToCheck]) - np.min(abs(Jet_eta[:jetsToCheck]))) / (np.max(abs(Jet_eta[:jetsToCheck])) - np.min(abs(Jet_eta[:jetsToCheck]))))
            scores = scores + (Jet_qgl[:jetsToCheck] - np.min(Jet_qgl[:jetsToCheck])) / (np.max(Jet_qgl[:jetsToCheck]) - np.min(Jet_qgl[:jetsToCheck]))
            scores =  scores + (Jet_btagDeepFlavB[:jetsToCheck] - np.min(Jet_btagDeepFlavB[:jetsToCheck])) / (np.max(Jet_btagDeepFlavB[:jetsToCheck]) - np.min(Jet_btagDeepFlavB[:jetsToCheck]))
            scores = scores + IsTrigJet
            taken = np.argsort(scores)[::-1][:4]


            
            if np.array_equal(np.sort(taken), np.arange(nJet)[m]):
                correct = correct + 1

    print(correct/matchedEvents)
    print(matchedEvents/totalEntries)
    return

if __name__ == "__main__":
    main()