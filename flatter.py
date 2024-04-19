import numpy as np
import pandas as pd
import glob
import uproot
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
def jetsSelector(jetsToCheck, Jet_muonIdx1, Jet_muonIdx2, Muon_isTriggering, Jet_pt, Jet_eta, Jet_qgl, Jet_btagDeepFlavB):
    # list of bools [nJet]
    IsTrigJet = []
    muonIdxs = []
    muonIdx = 999
    # list of scores based on pT (max), mass (max), |eta| (min), btag (max)
    scores = []
    
    for jetIdx in range(jetsToCheck):
        #jetTemp = ROOT.TLorentzVector(0.,0.,0.,0.)
        #jetTemp.SetPtEtaPhiM(Jet_pt[jetIdx], Jet_eta[jetIdx], Jet_phi[jetIdx], Jet_mass[jetIdx])
        #myRootJets.append(jetTemp)
        
        if Jet_muonIdx1[jetIdx]>-1:
            if (Muon_isTriggering[Jet_muonIdx1[jetIdx]]):
                IsTrigJet.append(100)
                muonIdxs.append(Jet_muonIdx1[jetIdx])
                continue
        if Jet_muonIdx2[jetIdx]>-1:
            if (Muon_isTriggering[Jet_muonIdx2[jetIdx]]):
                IsTrigJet.append(100)
                muonIdxs.append(Jet_muonIdx2[jetIdx])
            else:
                IsTrigJet.append(0)
        else:
            IsTrigJet.append(0)
    # temporary
    muonIdx = np.argmax(Muon_pt[muonIdxs])
    
    scores = (Jet_pt[:jetsToCheck] - np.min(Jet_pt[:jetsToCheck])) / (np.max(Jet_pt[:jetsToCheck]) - np.min(Jet_pt[:jetsToCheck]))
    scores = scores + (1- (abs(Jet_eta[:jetsToCheck]) - np.min(abs(Jet_eta[:jetsToCheck]))) / (np.max(abs(Jet_eta[:jetsToCheck])) - np.min(abs(Jet_eta[:jetsToCheck]))))
    scores = scores + (Jet_qgl[:jetsToCheck] - np.min(Jet_qgl[:jetsToCheck])) / (np.max(Jet_qgl[:jetsToCheck]) - np.min(Jet_qgl[:jetsToCheck]))
    scores =  scores + (Jet_btagDeepFlavB[:jetsToCheck] - np.min(Jet_btagDeepFlavB[:jetsToCheck])) / (np.max(Jet_btagDeepFlavB[:jetsToCheck]) - np.min(Jet_btagDeepFlavB[:jetsToCheck]))
    scores = scores + IsTrigJet
    taken = np.argsort(scores)[::-1][:4]

    return taken, muonIdx

def jetsBuilder(takenJets, Jet_pt, Jet_bReg2018, Jet_eta, Jet_phi, Jet_mass):
    myJets = []
    for jetIdx in takenJets:
        jetTemp = ROOT.TLorentzVector(0.,0.,0.,0.)
        jetTemp.SetPtEtaPhiM(Jet_pt[jetIdx]*Jet_bReg2018[jetIdx], Jet_eta[jetIdx], Jet_phi[jetIdx], Jet_mass[jetIdx])
        myJets.append(jetTemp)

    return myJets

def flatter():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:6]
    print("%d files to be used" %len(fileNames))
    totalEntries = 0
    HHmass = []
    
    for fileName in fileNames:
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries 
        totalEntries = totalEntries + maxEntries
        print("Entries : %d" %maxEntries)

        for ev in  range(maxEntries):
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
            #m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            #if np.sum(m)!=4:
            #    continue
            jetsToCheck = np.min((12, nJet))
            if nJet<4:
                continue
            takenJets, muonIdx = jetsSelector(jetsToCheck, Jet_muonIdx1, Jet_muonIdx2, Muon_isTriggering, Jet_pt, Jet_eta, Jet_qgl, Jet_btagDeepFlavB)
            myJets = jetsBuilder(takenJets, Jet_pt, Jet_bReg2018, Jet_eta, Jet_phi, Jet_mass)
            HHmass.append((myJets[0]+myJets[1]+myJets[2]+myJets[3]).M())
    
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/0000"
    fileNames = glob.glob(path+"/Data2018*.root")[:6]
    print("%d files to be used" %len(fileNames))
    totalEntries = 0
    HHmassBkg = []
    
    for fileName in fileNames:
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries 
        totalEntries = totalEntries + maxEntries
        print("Entries : %d" %maxEntries)

        for ev in  range(maxEntries):
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
            #m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            #if np.sum(m)!=4:
            #    continue
            jetsToCheck = np.min((12, nJet))
            if nJet<4:
                continue
            takenJets = jetsSelector(jetsToCheck, Jet_muonIdx1, Jet_muonIdx2, Muon_isTriggering, Jet_pt, Jet_eta, Jet_qgl, Jet_btagDeepFlavB)
            myJets = jetsBuilder(takenJets, Jet_pt, Jet_bReg2018, Jet_eta, Jet_phi, Jet_mass)
            HHmass.append((myJets[0]+myJets[1]+myJets[2]+myJets[3]).M())

    fig, ax = plt.subplots(1, 1)
    bins = np.linspace(100, 1200, 40)
    ax.hist(np.clip(HHmass, bins[0], bins[-1]),bins=bins, color='blue', histtype=u'step')
    ax.set_xlabel("m(HH) [GeV]")
    fig.savefig("/t3home/gcelotto/ggHH4b/HHmass.png", bbox_inches='tight')

    return

if __name__ == "__main__":
    flatter()