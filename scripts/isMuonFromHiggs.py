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
    # numerator and denominator for the ratio
    trigJetFromHiggs = 0
    trigJet = 0
    muonFromHiggs = 0

    for fileName in fileNames:
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries 
        totalEntries = totalEntries + maxEntries
        print("Entries : %d" %maxEntries)

        for ev in  range(maxEntries):
            GenPart_statusFlags         = branches["GenPart_statusFlags"][ev]
            GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
            GenPart_pdgId               = branches["GenPart_pdgId"][ev]

            Muon_genPartIdx             = branches["Muon_genPartIdx"][ev]
            GenPart_genPartIdxMother    = branches["GenPart_genPartIdxMother"][ev]
            GenPart_pdgId               = branches["GenPart_pdgId"][ev]

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
            def iterativeMother(GenPart_pdgId,  GenPart_genPartIdxMother, genPartIdx0):
                if genPartIdx0==-1:
                    return False
                if GenPart_genPartIdxMother[genPartIdx0] != -1:
                    if GenPart_pdgId[GenPart_genPartIdxMother[genPartIdx0]]==25:
                    
                        return True
                    else:
                        genPartIdx0 = GenPart_genPartIdxMother[genPartIdx0]
                    
                        result=iterativeMother(GenPart_pdgId,  GenPart_genPartIdxMother, genPartIdx0)
                        return result
                else:
                    return False
                
            def muonGenMatched(muon_genPartIdx):
                if muon_genPartIdx==-1:
                    return False
                else:
                    return True
                
            #def isMuonMatched(GenPart_pdgId,  GenPart_genPartIdxMother, genPartIdx0):


                
            matchedEvents=matchedEvents+1
            
            for jetIdx in range(nJet):
                # check if trig jets in matched jets
                trigJetFilled = False
                if Jet_muonIdx1[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx1[jetIdx]]):
                        # check the muon
                        muonFromHiggs = muonFromHiggs + iterativeMother(GenPart_pdgId,  GenPart_genPartIdxMother, Muon_genPartIdx[Jet_muonIdx1[jetIdx]])
                        #muonFromHiggs = muonFromHiggs + muonGenMatched(Muon_genPartIdx[Jet_muonIdx1[jetIdx]])
                        trigJet = trigJet + 1
                        trigJetFilled = True
                        if jetIdx in np.arange(nJet)[m]:
                            trigJetFromHiggs = trigJetFromHiggs + 1
                            continue
                if Jet_muonIdx2[jetIdx]>-1:
                    if (Muon_isTriggering[Jet_muonIdx2[jetIdx]]):
                        muonFromHiggs = muonFromHiggs + iterativeMother(GenPart_pdgId,  GenPart_genPartIdxMother, Muon_genPartIdx[Jet_muonIdx2[jetIdx]])
                        #muonFromHiggs = muonFromHiggs + muonGenMatched(Muon_genPartIdx[Jet_muonIdx2[jetIdx]])
                        # check the muon
                        if not trigJetFilled:
                            trigJet = trigJet + 1
                        if jetIdx in np.arange(nJet)[m]:
                            trigJetFromHiggs = trigJetFromHiggs + 1
                            continue
    print(trigJetFromHiggs/trigJet)
    print(muonFromHiggs, trigJet, muonFromHiggs/trigJet)

if __name__ == "__main__":
    main()
            


