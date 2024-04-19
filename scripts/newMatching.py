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
def main(nFiles):
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:nFiles]
    print("%d files to be used" %len(fileNames))
    matchedEvents, totalEntries = 0, 0, 
    completeMatching = 0
    partialMatching = 0
    partialNotWorking = 0

    data = {'mH1' : [],
            'mH2' : [],
            'lowestJetPt' : []}
    
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
                    j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
                    j2 = ROOT.TLorentzVector(0.,0.,0.,0.)
                    j3 = ROOT.TLorentzVector(0.,0.,0.,0.)
                    j4 = ROOT.TLorentzVector(0.,0.,0.,0.)

                    j1.SetPtEtaPhiM(Jet_pt[newM][0]*Jet_bReg2018[newM][0], Jet_eta[newM][0], Jet_phi[newM][0], Jet_mass[newM][0])
                    j2.SetPtEtaPhiM(Jet_pt[newM][1]*Jet_bReg2018[newM][1], Jet_eta[newM][1], Jet_phi[newM][1], Jet_mass[newM][1])
                    j3.SetPtEtaPhiM(Jet_pt[newM][2]*Jet_bReg2018[newM][2], Jet_eta[newM][2], Jet_phi[newM][2], Jet_mass[newM][2])
                    j4.SetPtEtaPhiM(Jet_pt[newM][3]*Jet_bReg2018[newM][3], Jet_eta[newM][3], Jet_phi[newM][3], Jet_mass[newM][3])
                    myJets = [j1, j2, j3, j4]

                    for i in np.arange(4):
                        for j in np.arange(4):
                            if i == j:
                                continue
                            elif GenJet_partonMotherIdx[Jet_genJetIdx[newM][i]] == GenJet_partonMotherIdx[Jet_genJetIdx[newM][j]]:
                                couple1 = [i, j]
                                continue
                    couple2=[]
                    for i in np.arange(4):
                        if i not in couple1:
                            couple2.append(i)
                    
                    #print(couple1, couple2)
                    if (myJets[couple1[0]] + myJets[couple1[1]]).Pt() < (myJets[couple2[0]] + myJets[couple2[1]]).Pt():
                        couple1, couple2 = couple2, couple1
                    assert (myJets[couple1[0]] + myJets[couple1[1]]).Pt() > (myJets[couple2[0]] + myJets[couple2[1]]).Pt(), "Assumption Failed in ev %d"
                    data['mH1'].append((myJets[couple1[0]] + myJets[couple1[1]]).M())
                    data['mH2'].append((myJets[couple2[0]] + myJets[couple2[1]]).M())
                    data['lowestJetPt'].append(np.min((myJets[0].Pt(), myJets[1].Pt(), myJets[2].Pt(), myJets[3].Pt())))
                    
                    partialMatching = partialMatching +1
                    m = newM
                elif np.sum(newM)!=4:
                    partialNotWorking = partialNotWorking + 1


    print("Partial Matching : ", partialMatching/totalEntries)
    print("Complete Matching : ", completeMatching/totalEntries)
    print("Partial not working : ", partialNotWorking/totalEntries)

    print("saving the df")
    df = pd.DataFrame(data)
    outPath = "/t3home/gcelotto/ggHH4b/flatData/df_mH1mH2_genMatched.parquet"
    df.to_parquet(outPath)

    return

if __name__ == "__main__":
    nFiles = int(sys.argv[1]) if len(sys.argv)>1 else None
    if nFiles == None:
        sys.exit("Specify nFiles")
    main(nFiles=nFiles)