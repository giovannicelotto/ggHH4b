import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:5]
    print("%d files to be used" %len(fileNames))
    totalEntries = 0
    numberOfMatched=[]
    numberOfGen=[]

    fourthJetPt, fourthJetEtaPhi = [], []
    
    for fileName in fileNames:
        print("Opening ", fileName)
        f = uproot.open(fileName)
        tree = f['Events']
        branches = tree.arrays()
        maxEntries = tree.num_entries 
        totalEntries = totalEntries + maxEntries
        print("Entries : %d" %maxEntries)

        for ev in  range(maxEntries):
            GenPart_pt                  = branches["GenPart_pt"][ev]
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
            
            mReco = (Jet_genJetIdx>-1) & (abs(GenJet_partonFlavour[Jet_genJetIdx])==5) & (GenJet_partonMotherPdgId[Jet_genJetIdx]==25)
            mGen =  (abs(GenJet_partonFlavour)==5) & (GenJet_partonMotherPdgId==25)
            numberOfMatched.append(np.sum(mReco))
            numberOfGen.append(np.sum(mGen))

            if np.sum(mReco)!=3:
                continue

            input("%d Quarks"%ev)
            mQuark = (GenPart_pdgId[GenPart_genPartIdxMother]==25) & ((GenPart_statusFlags[GenPart_genPartIdxMother] & 8192)!=0)
            print(GenPart_pt[mQuark], GenPart_eta[mQuark])
            print(GenJet_pt[mGen], GenJet_eta[mGen])
            

            flags = ((GenPart_statusFlags & 8192)!=0) & (GenPart_pdgId==25)
            h1 = ROOT.TLorentzVector(0.,0.,0.,0.)
            h2 = ROOT.TLorentzVector(0.,0.,0.,0.)
            j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
            j2 = ROOT.TLorentzVector(0.,0.,0.,0.)
            j3 = ROOT.TLorentzVector(0.,0.,0.,0.)

            h1.SetPtEtaPhiM(GenPart_pt[flags][0], GenPart_eta[flags][0], GenPart_phi[flags][0], GenPart_mass[flags][0])
            h2.SetPtEtaPhiM(GenPart_pt[flags][1], GenPart_eta[flags][1], GenPart_phi[flags][1], GenPart_mass[flags][1])
            j1.SetPtEtaPhiM(GenJet_pt[mGen][0], GenJet_eta[mGen][0], GenJet_phi[mGen][0], GenJet_mass[mGen][0])
            j2.SetPtEtaPhiM(GenJet_pt[mGen][1], GenJet_eta[mGen][1], GenJet_phi[mGen][1], GenJet_mass[mGen][1])
            j3.SetPtEtaPhiM(GenJet_pt[mGen][2], GenJet_eta[mGen][2], GenJet_phi[mGen][2], GenJet_mass[mGen][2])

            fourthJetPt.append((h1+h2-j1-j2-j3).Pt())
            fourthJetEtaPhi.append(((h1+h2-j1-j2-j3).Eta(), ((h1+h2-j1-j2-j3).Phi())))

    fig, ax = plt.subplots(1, 2, figsize=(9, 6))
    bins=np.linspace(0, 100, 6)
    c = np.histogram(fourthJetPt, bins=bins)[0]
    c = c/np.sum(c)
    ax[0].hist(bins[:-1], bins=bins, weights=c, color='blue', histtype=u'step')
    ax[0].set_xlabel("Jet pt [GeV]")

    x_values, y_values = zip(*fourthJetEtaPhi)
    c = ax[1].hist2d(x_values, y_values, bins=(10, 10), cmap=plt.cm.viridis)[0]
    ax[1].set_xlabel("Jet eta")
    ax[1].set_ylabel("Jet phi")

    outName = "/t3home/gcelotto/ggHH4b/fourthJetInfo.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in", outName)


    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    bins=np.linspace(0, 5, 6)
    c = np.histogram(numberOfMatched, bins=bins)[0]
    c = c/np.sum(c)
    ax.hist(bins[:-1], bins=bins-0.5, weights=c, color='blue', histtype=u'step')
    ax.set_xlabel("N Jet Matched")
    outName = "/t3home/gcelotto/ggHH4b/numMatched.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in", outName)

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    bins=np.linspace(0, 5, 6)
    c = np.histogram(numberOfGen, bins=bins)[0]
    c = c/np.sum(c)
    ax.hist(bins[:-1], bins=bins-0.5, weights=c, color='blue', histtype=u'step')
    ax.set_xlabel("N Gen Jet From Higgs")
    outName = "/t3home/gcelotto/ggHH4b/numHiggsDaughtersGen.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in", outName)
    return

if __name__ == "__main__":
    main()