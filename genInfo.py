import numpy as np
import uproot
import pandas as pd
import glob
import sys
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
def main(saveAndLoad, nFiles):
    if saveAndLoad:
        path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
        fileNames = glob.glob(path+"/ggHH4b*.root")[:nFiles]
        print("%d files to be used" %len(fileNames))
        matchedEvents, totalEntries = 0, 0
        features_=[]
        bkgFeatures_ = []
        for fileName in fileNames:
            f = uproot.open(fileName)
            tree = f['Events']
            branches = tree.arrays()
            maxEntries = tree.num_entries 
            totalEntries = totalEntries + maxEntries
            print("Entries : %d" %maxEntries)

            for ev in  range(maxEntries):
                eventFeatures_ =  {}
                bkgEventJetsFeatures =  {}

                if (ev%(int(maxEntries/20))==0):
                    sys.stdout.write('\r')
                    sys.stdout.write("%d%%"%(ev/maxEntries*100))
                    sys.stdout.flush()


                GenJet_partonFlavour        = branches["GenJet_partonFlavour"][ev]
                GenJet_partonMotherIdx      = branches["GenJet_partonMotherIdx"][ev]
                GenJet_partonMotherPdgId    = branches["GenJet_partonMotherPdgId"][ev]
            # Reco Jets
                nJet                        = branches["nJet"][ev]
                Jet_eta                     = branches["Jet_eta"][ev]
                Jet_pt                      = branches["Jet_pt"][ev]
                Jet_phi                     = branches["Jet_phi"][ev]
                Jet_mass                    = branches["Jet_mass"][ev]
                Jet_bReg2018                = branches["Jet_bReg2018"][ev]
                Jet_genJetIdx               = branches["Jet_genJetIdx"][ev]
                Jet_btagDeepFlavB           = branches["Jet_btagDeepFlavB"][ev]
                Jet_muonIdx1                = branches["Jet_muonIdx1"][ev]
                Jet_muonIdx2                = branches["Jet_muonIdx2"][ev]
                Muon_isTriggering           = branches["Muon_isTriggering"][ev]
                GenJet_pt                   = branches["GenJet_pt"][ev]
                GenJet_eta                  = branches["GenJet_eta"][ev]
                GenJet_phi                  = branches["GenJet_phi"][ev]
                GenJet_mass                 = branches["GenJet_mass"][ev]

# mask events with 4 reco mathced jets, from Higgs, b-jets
                m = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
                if np.sum(m)!=4:
                    continue
                matchedEvents=matchedEvents+1

# h1 is the higgs with the leading jet
                h1_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]==GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
                h2_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]!=GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]



                h1j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
                h1j2 = ROOT.TLorentzVector(0.,0.,0.,0.)
                h2j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
                h2j2 = ROOT.TLorentzVector(0.,0.,0.,0.)

                h1j1.SetPtEtaPhiM(GenJet_pt[h1_GenJetIdx[0]], GenJet_eta[h1_GenJetIdx[0]], GenJet_phi[h1_GenJetIdx[0]], GenJet_mass[h1_GenJetIdx[0]])
                h1j2.SetPtEtaPhiM(GenJet_pt[h1_GenJetIdx[1]], GenJet_eta[h1_GenJetIdx[1]], GenJet_phi[h1_GenJetIdx[1]], GenJet_mass[h1_GenJetIdx[1]])
                h2j1.SetPtEtaPhiM(GenJet_pt[h2_GenJetIdx[0]], GenJet_eta[h2_GenJetIdx[0]], GenJet_phi[h2_GenJetIdx[0]], GenJet_mass[h2_GenJetIdx[0]])
                h2j2.SetPtEtaPhiM(GenJet_pt[h2_GenJetIdx[1]], GenJet_eta[h2_GenJetIdx[1]], GenJet_phi[h2_GenJetIdx[1]], GenJet_mass[h2_GenJetIdx[1]])

                Jet_muonIdx1


                eventFeatures_['h1j1pt'] = (h1j1.Pt())
                eventFeatures_['h1j1eta'] = (h1j1.Eta())
                eventFeatures_['h1j1mass'] = (h1j1.M())
                eventFeatures_['h1j2pt'] = (h1j2.Pt())
                eventFeatures_['h1j2eta'] = (h1j2.Eta())
                eventFeatures_['h1j2mass'] = (h1j2.M())
                eventFeatures_['h2j1pt'] = (h2j1.Pt())
                eventFeatures_['h2j1eta'] = (h2j1.Eta())
                eventFeatures_['h2j1mass'] = (h2j1.M())
                eventFeatures_['h2j2pt'] = (h2j2.Pt())
                eventFeatures_['h2j2eta'] = (h2j2.Eta())
                eventFeatures_['h2j2mass'] = (h2j2.M())
                eventFeatures_['h2j2pt'] = (h2j2.Pt())
                eventFeatures_['h2j2eta'] = (h2j2.Eta())
                if (h1j1+h1j2).Pt() > (h2j1+h2j2).Pt():
                    eventFeatures_['h1_dR'] = (h1j1.DeltaR(h1j2))
                    eventFeatures_['h2_dR'] = (h2j1.DeltaR(h2j2))
  
                else:
                    eventFeatures_['h1_dR'] = (h2j1.DeltaR(h2j2))
                    eventFeatures_['h2_dR'] = (h1j1.DeltaR(h1j2))
                    
                
                bkgJetsIdx = np.arange(nJet)[~ m]
                for idx in bkgJetsIdx:
                    bkgEventJetsFeatures = {}
                    bkgEventJetsFeatures['jetpt'] = Jet_pt[idx]
                    bkgEventJetsFeatures['jeteta'] = Jet_eta[idx]
                    bkgEventJetsFeatures['jetmass'] = Jet_mass[idx]
                    bkgFeatures_.append(bkgEventJetsFeatures)
                features_.append(eventFeatures_)
        features = pd.DataFrame(features_)
        bkgFeatures = pd.DataFrame(bkgFeatures_)
        
        del features_, bkgFeatures_
        features.to_csv("/t3home/gcelotto/ggHH4b/genFeatures.csv")
        bkgFeatures.to_csv("/t3home/gcelotto/ggHH4b/bkgGenFeatures.csv")
    print("loading the df")
    features    =   pd.read_csv("/t3home/gcelotto/ggHH4b/genFeatures.csv")
    bkgFeatures =   pd.read_csv("/t3home/gcelotto/ggHH4b/bkgGenFeatures.csv")
    features=features.iloc[:,1:]
    bkgFeatures=bkgFeatures.iloc[:,1:]

    #print(bkgFeatures.shape)
    #print(bkgFeatures.describe())
    print(bkgFeatures.head(3))
    
    #
    #print(bkgFeatures.head())
    nRow, nCol = 7, 3
    fig, ax = plt.subplots(nRow, nCol, constrained_layout=True)
    for i in range(nRow):
        for j in range(nCol):
            if i*nCol+j>=features.shape[1]:
                break
            c, b_ = np.histogram(features.iloc[:,i*nCol+j], bins=20)[:2]
            ax[i, j].hist(b_[:-1], bins=b_, weights=c)
    outName = "/t3home/gcelotto/ggHH4b/hh_genFeatures.png"
    fig.savefig(outName, bbox_inches='tight')
    print("saving in ", outName)

    nRow, nCol = 2, 3
    fig, ax = plt.subplots(nRow, nCol, constrained_layout=True, figsize=(10, 6))
    
    bins_pt     = np.linspace(0, 200, 30)
    bins_eta    = np.linspace(-5, 5, 30)
    bins_mass   = np.linspace(0, 20, 30)
    counts_pt   = np.histogram(features.iloc[:,[0, 3, 6, 9]].values.flatten(), bins=bins_pt)[0]
    counts_eta  = np.histogram(features.iloc[:,[1, 4, 7, 10]].values.flatten(), bins=bins_eta)[0]
    counts_mass = np.histogram(features.iloc[:,[2, 5, 8, 11]].values.flatten(), bins=bins_mass)[0]
    counts_pt   = counts_pt/np.sum(counts_pt)
    counts_eta  = counts_eta/np.sum(counts_eta)
    counts_mass = counts_mass/np.sum(counts_mass)
    ax[0, 0].hist(bins_pt[:-1], bins=bins_pt, weights=counts_pt, histtype=u'step', color='blue', label='Signal')
    ax[0, 1].hist(bins_eta[:-1], bins=bins_eta, weights=counts_eta, histtype=u'step', color='blue', label='Signal')
    ax[0, 2].hist(bins_mass[:-1], bins=bins_mass, weights=counts_mass, histtype=u'step', color='blue', label='Signal')

    counts_pt   = np.histogram(np.clip(bkgFeatures.jetpt, bins_pt[0], bins_pt[-1]), bins=bins_pt)[0]
    counts_eta  = np.histogram(np.clip(bkgFeatures.jeteta, bins_eta[0], bins_eta[-1]), bins=bins_eta)[0]
    counts_mass = np.histogram(np.clip(bkgFeatures.jetmass, bins_mass[0], bins_mass[-1]), bins=bins_mass)[0]
    counts_pt   = counts_pt/np.sum(counts_pt)
    counts_eta  = counts_eta/np.sum(counts_eta)
    counts_mass = counts_mass/np.sum(counts_mass)
    ax[0, 0].hist(bins_pt[:-1], bins=bins_pt, weights=counts_pt, histtype=u'step', color='red', label='Other Jets')
    ax[0, 1].hist(bins_eta[:-1], bins=bins_eta, weights=counts_eta, histtype=u'step', color='red', label='Other Jets')
    ax[0, 2].hist(bins_mass[:-1], bins=bins_mass, weights=counts_mass, histtype=u'step', color='red', label='Other Jets')
    ax[0, 0].set_yscale('log')

    
    
    bins_dR      = np.linspace(0, 5, 20)
    counts_dR1   = np.histogram(features.h1_dR, bins=bins_dR)[0]
    counts_dR2   = np.histogram(features.h2_dR, bins=bins_dR)[0]
    counts_dR1   = counts_dR1/np.sum(counts_dR1)
    counts_dR2   = counts_dR2/np.sum(counts_dR2)
    ax[1, 0].hist(bins_dR[:-1], bins=bins_dR, weights=counts_dR1, histtype=u'step', color='blue', label='Signal')
    ax[1, 1].hist(bins_dR[:-1], bins=bins_dR, weights=counts_dR2, histtype=u'step', color='blue', label='Signal')


    
    outName = "/t3home/gcelotto/ggHH4b/hh_genInfo.png"
    fig.savefig(outName, bbox_inches='tight')
    print("saving in ", outName)

    return

if __name__ == "__main__":
    if len(sys.argv)==2:
        saveAndLoad=bool(int(sys.argv[1]))
        main(saveAndLoad, nFiles=None)
    if len(sys.argv)>2:
        saveAndLoad=bool(int(sys.argv[1]))
        nFiles=int(sys.argv[2])
        main(saveAndLoad, nFiles)