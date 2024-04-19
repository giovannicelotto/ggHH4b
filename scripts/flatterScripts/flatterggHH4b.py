import numpy as np
import uproot
import glob
import awkward as ak
import pandas as pd
import sys
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import re
import os
hep.style.use("CMS")

def jetsSelector(jetsToCheck, Muon_isTriggering, Jet_muonIdx1, Jet_muonIdx2, Jet_pt, Jet_eta, Jet_qgl, Jet_btagDeepFlavB, Muon_pt):

    Jet_qgl_minBound = ak.where(Jet_qgl<0, 0, Jet_qgl)
    IsTrigJet, muonIdxs=[], []
    for jetIdx in range(jetsToCheck):
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
    #take the muonIdx of the leading muon
    muonIdx = np.argmax(Muon_pt[muonIdxs])

    p0, p1, p2, p3 = 0.3034,    0.1752,  0.1849,     0.6697 
    scores = p0*(Jet_pt[:jetsToCheck] - np.min(Jet_pt[:jetsToCheck])) / (np.max(Jet_pt[:jetsToCheck]) - np.min(Jet_pt[:jetsToCheck]))
    scores = scores + p1*(1- (abs(Jet_eta[:jetsToCheck]) - np.min(abs(Jet_eta[:jetsToCheck]))) / (np.max(abs(Jet_eta[:jetsToCheck])) - np.min(abs(Jet_eta[:jetsToCheck]))))
    scores = scores + p2*((Jet_qgl_minBound[:jetsToCheck] - np.min(Jet_qgl_minBound[:jetsToCheck])) / (np.max(Jet_qgl_minBound[:jetsToCheck]) - np.min(Jet_qgl_minBound[:jetsToCheck])+1e-10))
    scores = scores + p3*((Jet_btagDeepFlavB[:jetsToCheck] - np.min(Jet_btagDeepFlavB[:jetsToCheck])) / (np.max(Jet_btagDeepFlavB[:jetsToCheck]) - np.min(Jet_btagDeepFlavB[:jetsToCheck])))
    scores = scores + IsTrigJet
    if np.max(scores)<100: #no trig jets
        return -999, -999
    else:
        taken = list(np.argsort(scores)[::-1][:4])
        return taken, muonIdx

def getRelativePhi(v1, v2):
    dPhi = v1.Phi() - v2.Phi()
    dPhi = abs(dPhi - 2*np.pi*(dPhi > np.pi) + 2*np.pi*(dPhi< -np.pi))
    return np.float32(dPhi)

def pairRelations(v1, v2, noPhi=False):
            dPhi = v1.Phi() - v2.Phi()
            dPhi = abs(dPhi - 2*np.pi*(dPhi > np.pi) + 2*np.pi*(dPhi< -np.pi))

            dEta = abs(v1.Eta() - v2.Eta())

            dR = v1.DeltaR(v2)
            tau = np.arctan(abs(dPhi)/(dEta + 0.0000001))
            if noPhi:
                return [np.float32(dEta), np.float32(dR), np.float32(tau)]
            return [np.float32(dPhi), np.float32(dEta), np.float32(dR), np.float32(tau)]

def distanceHiggsPlane(m1, m2):
                x0 = (m1 + m2)/2
                y0 = x0
                dist = np.sqrt((m1-x0)**2+(m2-y0)**2)
                return dist


def treeFlatten(fileName, maxEntries, isMC):
    f = uproot.open(fileName)
    tree = f['Events']
    branches = tree.arrays()
    maxEntries = tree.num_entries if maxEntries==-1 else maxEntries
    print("Entries : %d"%(maxEntries))
    file_ =[]
    histPath = "/t3home/gcelotto/ggHbb/trgMu_scale_factors.root"
    f = ROOT.TFile(histPath, "READ")
    hist = f.Get("hist_scale_factor")
    for ev in  range(maxEntries):
        
        features_ = []
        if (ev%(int(maxEntries/20))==0):
            sys.stdout.write('\r')
            # the exact output you're looking for:
            sys.stdout.write("%d%%"%(ev/maxEntries*100))
            sys.stdout.flush()
    # Reco Jets
        nJet                        = branches["nJet"][ev]
        Jet_pt                      = branches["Jet_pt"][ev]
        Jet_eta                     = branches["Jet_eta"][ev]
        Jet_phi                     = branches["Jet_phi"][ev]
        Jet_mass                    = branches["Jet_mass"][ev]
        Jet_bReg2018                = branches["Jet_bReg2018"][ev]
        Jet_nConstituents           = branches["Jet_nConstituents"][ev]
        
        Jet_btagDeepFlavB           = branches["Jet_btagDeepFlavB"][ev]
        
        Jet_muonIdx1                = branches["Jet_muonIdx1"][ev]
        Jet_muonIdx2                = branches["Jet_muonIdx2"][ev]
        Muon_isTriggering           = branches["Muon_isTriggering"][ev]
        Jet_qgl                     = branches["Jet_qgl"][ev]
        Jet_nMuons                  = branches["Jet_nMuons"][ev]
        Jet_nElectrons              = branches["Jet_nElectrons"][ev]

        Muon_pt                     = branches["Muon_pt"][ev]
        Muon_eta                    = branches["Muon_eta"][ev]
        Muon_phi                    = branches["Muon_phi"][ev]
        Muon_mass                   = branches["Muon_mass"][ev]

        
# STEP 1 
# choose 4 jets
        if nJet<4:
            continue
        
        jetsToCheck = np.min((12, nJet))
        muonIdx = -999
        taken = -999
        taken, muonIdx = jetsSelector(jetsToCheck, Muon_isTriggering, Jet_muonIdx1, Jet_muonIdx2, Jet_pt, Jet_eta, Jet_qgl, Jet_btagDeepFlavB, Muon_pt)
        
        if (muonIdx==-999) & (taken==-999):
            #jets selector return -999 when there is no triggering muon
            continue
        # taken are sorted by scores. j1 will always have the triggering muon
        

        j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
        j2 = ROOT.TLorentzVector(0.,0.,0.,0.)
        j3 = ROOT.TLorentzVector(0.,0.,0.,0.)
        j4 = ROOT.TLorentzVector(0.,0.,0.,0.)
        h1 = ROOT.TLorentzVector(0.,0.,0.,0.)
        h2 = ROOT.TLorentzVector(0.,0.,0.,0.)

        j1.SetPtEtaPhiM(Jet_pt[taken[0]]*Jet_bReg2018[taken[0]], Jet_eta[taken[0]], Jet_phi[taken[0]], Jet_mass[taken[0]])
        j2.SetPtEtaPhiM(Jet_pt[taken[1]]*Jet_bReg2018[taken[1]], Jet_eta[taken[1]], Jet_phi[taken[1]], Jet_mass[taken[1]])
        j3.SetPtEtaPhiM(Jet_pt[taken[2]]*Jet_bReg2018[taken[2]], Jet_eta[taken[2]], Jet_phi[taken[2]], Jet_mass[taken[2]])
        j4.SetPtEtaPhiM(Jet_pt[taken[3]]*Jet_bReg2018[taken[3]], Jet_eta[taken[3]], Jet_phi[taken[3]], Jet_mass[taken[3]])
        
            
            #print("Muon is triggering ", Muon_isTriggering[Jet_muonIdx1[taken[0]]])
            #print(taken[0], muonIdx)
            #print(ev)
        
        # 1-2 3-4
        
        dist1 = distanceHiggsPlane(m1=(j1+j2).M(), m2=(j3+j4).M())
        dist2 = distanceHiggsPlane(m1=(j1+j3).M(), m2=(j2+j4).M())
        dist3 = distanceHiggsPlane(m1=(j1+j4).M(), m2=(j2+j3).M())
        m1, m2 = np.array([(j1+j2).M(), (j1+j3).M(), (j1+j4).M()]), np.array([(j3+j4).M(), (j2+j4).M(), (j2+j3).M()])
        minDistance = np.min((dist1, dist2, dist3))
        minIndex    = np.argmin((dist1, dist2, dist3))
        if abs(minDistance - sorted([dist1, dist2, dist3])[1]) < 25:
            minIndex    = np.argmin((m1+m2)/2)
            minDistance = [dist1, dist2, dist3][np.argmin((m1+m2)/2)]
        
        
        if minIndex   == 0: # j1 and j2
            pass
        elif minIndex == 1: # j1 and j3
            j2, j3 = j3, j2
            taken[1], taken[2] = taken[2], taken[1]
        elif minIndex == 2: # j1 and j4
            j2, j4 = j4, j2
            taken[1], taken[3] = taken[3], taken[1]
        # j1 is now always paired to j2
        # j3 and j4 are automatically matched
        # order j3 and j4 by pt
        
        # taken_0 = j1
        # taken_1 = j2
        # taken_2 = j3
        # taken_3 = j4
            
        if j3.Pt()<j4.Pt():
            j3, j4 = j4, j3
        
        # now order the higgs h1 is the leading
        if (j1+j2).Pt()>(j3+j4).Pt():
            h1 = j1 + j2
            h2 = j3 + j4
        else:
            h1 = j3 + j4
            h2 = j1 + j2
        features_.append((h1+h2).Pt())   
        features_.append((h1+h2).Eta())   
        features_.append(getRelativePhi((h1+h2), h1))   
        features_.append((h1+h2).M())   
        features_.append(h1.Pt())
        features_.append(h1.Eta())
        features_.append(h1.Phi())
        features_.append(h1.M())

        features_.append(h2.Pt())
        features_.append(h2.Eta())
        features_.append(getRelativePhi(h1, h2))
        features_.append(h2.M())

        features_.append(j1.Pt())
        features_.append(j1.Eta())
        features_.append(getRelativePhi(h1, j1))
        features_.append(j1.M())

        features_.append(j2.Pt())
        features_.append(j2.Eta())
        features_.append(getRelativePhi(h1, j2))
        features_.append(j2.M())
        

        features_.append(j3.Pt())
        features_.append(j3.Eta())
        features_.append(getRelativePhi(h1, j3))
        features_.append(j3.M())

        features_.append(j4.Pt())
        features_.append(j4.Eta())
        features_.append(getRelativePhi(h1, j4))
        features_.append(j4.M())

        features_.append(Jet_btagDeepFlavB[taken[0]])
        features_.append(Jet_btagDeepFlavB[taken[1]])
        features_.append(Jet_btagDeepFlavB[taken[2]])
        features_.append(Jet_btagDeepFlavB[taken[3]])

        features_.append(Jet_qgl[taken[0]])
        features_.append(Jet_qgl[taken[1]])
        features_.append(Jet_qgl[taken[2]])
        features_.append(Jet_qgl[taken[3]])

        features_.append(Jet_nConstituents[taken[0]])
        features_.append(Jet_nConstituents[taken[1]])
        features_.append(Jet_nConstituents[taken[2]])
        features_.append(Jet_nConstituents[taken[3]])

        features_.append(Jet_nMuons[taken[0]])
        features_.append(Jet_nMuons[taken[1]])
        features_.append(Jet_nMuons[taken[2]])
        features_.append(Jet_nMuons[taken[3]])

        features_.append(Jet_nElectrons[taken[0]])
        features_.append(Jet_nElectrons[taken[1]])
        features_.append(Jet_nElectrons[taken[2]])
        features_.append(Jet_nElectrons[taken[3]])

        features_ = features_ + pairRelations(h1, h2, noPhi=True)
        features_ = features_ + pairRelations(j1, j2)
        features_ = features_ + pairRelations(j3, j4)

        features_ = features_ + pairRelations(h1, j1, noPhi=True)
        features_ = features_ + pairRelations(h1, j2, noPhi=True)
        features_ = features_ + pairRelations(h2, j3)
        features_ = features_ + pairRelations(h2, j4)

        features_ = features_ + pairRelations(j1, j3)
        features_ = features_ + pairRelations(j1, j4)
        features_ = features_ + pairRelations(j2, j3)
        features_ = features_ + pairRelations(j2, j4)

        ht=np.sum(Jet_pt*Jet_bReg2018)
        features_.append(ht)
        myMuon = ROOT.TLorentzVector(0.,0.,0.,0.)
        
        myMuon.SetPtEtaPhiM(Muon_pt[muonIdx], Muon_eta[muonIdx], Muon_phi[muonIdx], Muon_mass[muonIdx])
        features_.append(myMuon.Pt())
        features_.append(myMuon.Eta())
        features_.append(getRelativePhi(myMuon, h1))
        file_.append(features_)
    return file_
    

def saveData(isMC, nFiles, maxEntries):
    pathSignal      = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    pathBkg         = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/000*"
    outFolderBkg    = "/t3home/gcelotto/bbar_analysis/flatData/selectedCandidates/data"
    outFolderSignal = "/t3home/gcelotto/bbar_analysis/flatData/selectedCandidates/ggHTrue"
    
    T3FolderSignal  = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData"
    T3FolderBkg     = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/flatDataForggHH4b"
    
    
    path = pathSignal if isMC else pathBkg
    
    
    fileNames = glob.glob(path+"/*.root")
    if nFiles > len(fileNames):
        print("nFiles > len(fileNames)\nSetting nFiles = len(fileNames) = %d" %len(fileNames))
        nFiles = len(fileNames)
    elif nFiles == -1:
        print("nFiles = -1\nSetting nFiles = len(fileNames) = %d" %len(fileNames))
        nFiles = len(fileNames)
    print("%d files to be used"%len(fileNames)) 

    doneFiles = 0
    for fileName in fileNames:
        if doneFiles == nFiles:
            sys.exit("Reached %d files\nExit"%nFiles)   
        outFolder = outFolderSignal if isMC else outFolderBkg
        T3Folder = T3FolderSignal   if isMC else T3FolderBkg
        fileNumber = re.search(r'\D(\d{1,4})\.\w+$', fileName).group(1)

        outName = "/ggHH4b_%s.parquet"%(fileNumber) if isMC else "/DataForggHH4b_%s.parquet"%(fileNumber)
        if os.path.exists(T3Folder +outName):
            # if you already saved this file skip
            print("File %s already present in T3\n" %fileNumber)
            continue
        if os.path.exists(outFolder + outName):
            # if you already saved this file skip
            print("File %s already present in bbar_analysis\n" %fileNumber)
            continue
        
        print("\nOpening ", (fileNames.index(fileName)+1), "/", len(fileNames), " path:", fileName, "...")
        df = pd.DataFrame()
        try:
            df.to_csv(T3Folder+outName, index=False, header=False)
            
        except:
            T3Folder="/scratch"
            df.to_csv(outFolder+outName, index=False, header=False)
        fileData = treeFlatten(fileName=fileName, maxEntries=maxEntries, isMC=isMC)
        print("Saving entries %d" %len(fileData))
        if T3Folder=="/scratch":
            print("Trying to remove ", outFolder+outName)
            os.remove(outFolder+outName)
            with open("/t3home/gcelotto/slurm/output/files.txt", 'a') as file:
                file.write(outName[1:]+"\n")
                print("Written %s in the file"%outName[1:])
        
        df=pd.DataFrame(fileData)
        featureNames = ['hh_pt', 'hh_eta', 'hh_RPhi','hh_mass',
                        'h1_pt', 'h1_eta', 'h1_phi', 'h1_mass',
                        'h2_pt', 'h2_eta', 'h2_Rphi', 'h2_mass',
                        'j1_pt', 'j1_eta', 'j1_Rphi', 'j1_mass',
                        'j2_pt', 'j2_eta', 'j2_Rphi', 'j2_mass',
                        'j3_pt', 'j3_eta', 'j3_Rphi', 'j3_mass',
                        'j4_pt', 'j4_eta', 'j4_Rphi', 'j4_mass',
                        'j1_btag',          'j2_btag',          'j3_btag',          'j4_btag',
                        'j1_qgl',           'j2_qgl',           'j3_qgl',           'j4_qgl',
                        'j1_nConstituents', 'j2_nConstituents', 'j3_nConstituents', 'j4_nConstituents',
                        'j1_nMuons',        'j2_nMuons',        'j3_nMuons',        'j4_nMuons',
                        'j1_nElectrons',    'j2_nElectrons',    'j3_nElectrons',    'j4_nElectrons',
                                     'dEta_h1h2', 'dR_h1h2', 'tau_h1h2',
                        'dPhi_j1j2', 'dEta_j1j2', 'dR_j1j2', 'tau_j1j2',
                        'dPhi_j3j4', 'dEta_j3j4', 'dR_j3j4', 'tau_j3j4',
                                     'dEta_h1j1', 'dR_h1j1', 'tau_h1j1',
                                     'dEta_h1j2', 'dR_h1j2', 'tau_h1j2',
                        'dPhi_h2j3', 'dEta_h2j3', 'dR_h2j3', 'tau_h2j3',
                        'dPhi_h2j4', 'dEta_h2j4', 'dR_h2j4', 'tau_h2j4',
                        
                        'dPhi_j1j3', 'dEta_j1j3', 'dR_j1j3', 'tau_j1j3',
                        'dPhi_j1j4', 'dEta_j1j4', 'dR_j1j4', 'tau_j1j4',
                        'dPhi_j2j3', 'dEta_j2j3', 'dR_j2j3', 'tau_j2j3',
                        'dPhi_j2j4', 'dEta_j2j4', 'dR_j2j4', 'tau_j2j4',
                        'ht', 'muon_pt', 'muon_eta', 'muon_Rphi']
        df.columns = featureNames
        df.to_parquet(T3Folder + outName)
        doneFiles = doneFiles +1
        #np.save(, np.array(fileData, np.float32))
        print("Saved in T3", T3Folder + outName)
    return 
     

if __name__=="__main__":
    isMC        = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    nFiles      = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    maxEntries  = int(sys.argv[3]) if len(sys.argv) > 3 else -1
    saveData(isMC, nFiles, maxEntries)
