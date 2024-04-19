import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep
import ROOT
hep.style.use("CMS")
def main():
    path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
    fileNames = glob.glob(path+"/ggHH4b*.root")[:10]
    print("%d files to be used" %len(fileNames))
    totalEntries = 0
    genHiggsPt, angDistance, etas = [], [], []
    
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
            Jet_bReg2018                = branches["Jet_bReg2018"][ev]
            Jet_genJetIdx               = branches["Jet_genJetIdx"][ev]
            Jet_btagDeepFlavB           = branches["Jet_btagDeepFlavB"][ev]
            GenJet_pt                   = branches["GenJet_pt"][ev]
            GenJet_eta                  = branches["GenJet_eta"][ev]
            GenJet_phi                  = branches["GenJet_phi"][ev]
            GenJet_mass                 = branches["GenJet_mass"][ev]
            
            flags = ((GenPart_statusFlags & 8192)!=0) & (GenPart_pdgId==25)  # higgs, isLastCopy()=True
            assert np.sum(flags)==2, "Number of Higgs is not 2"
            h1 = ROOT.TLorentzVector(0.,0.,0.,0.)
            h2 = ROOT.TLorentzVector(0.,0.,0.,0.)

            h1.SetPtEtaPhiM(GenPart_pt[flags][0], GenPart_eta[flags][0], GenPart_phi[flags][0], GenPart_mass[flags][0])
            h2.SetPtEtaPhiM(GenPart_pt[flags][1], GenPart_eta[flags][1], GenPart_phi[flags][1], GenPart_mass[flags][1])
            
            # swap so that h1 is leading
            if h1.Pt()<h2.Pt():
                h1, h2= h2, h1

            
            
            genHiggsPt.append((h1.Pt(), h2.Pt()))

            deltaEta = abs(h1.Eta() - h2.Eta())
            deltaPhi = abs(h1.DeltaPhi(h2))
            angDistance.append((deltaEta, deltaPhi))

            etas.append((h1.Eta(), h2.Eta()))
            

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    binsPt = np.linspace(0, 900, 40)
    pt1, pt2 = zip(*genHiggsPt)
    hist, xedges, yedges = np.histogram2d(pt1, pt2, bins=[binsPt, binsPt])
    im = ax.imshow(hist.T, cmap=plt.cm.viridis, norm=LogNorm(), origin='lower', extent=(binsPt.min(), binsPt.max(), binsPt.min(), binsPt.max()))
    cbar = fig.colorbar(im, ax=ax, label='Events', shrink=0.8)
    ax.set_xlabel("Leading Higgs Pt [GeV]")
    ax.set_ylabel("Subleading Higgs Pt [GeV]")
# diagonal line
    x_values = [0, binsPt[-1]]
    y_values = [0, binsPt[-1]]
    ax.plot(x_values, y_values, color='red', linestyle='--', label='')
    outName = "/t3home/gcelotto/ggHH4b/plots/genHiggsPt.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in", outName)

# NEW PLOT delta Eta Phi
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    binsEta = np.linspace(0, 5, 30)
    binsPhi = np.linspace(0, np.pi, 30)
    deltaEta, deltaPhi = zip(*angDistance)
    hist, xedges, yedges = np.histogram2d(deltaEta, deltaPhi, bins=[binsEta, binsPhi])
    im = ax.imshow(hist.T, cmap=plt.cm.viridis, norm=LogNorm(), origin='lower', aspect='auto', extent=(binsEta.min(), binsEta.max(), binsPhi.min(), binsPhi.max()))
    cbar = fig.colorbar(im, ax=ax, label='Events', shrink=1)  # Add color bar
    ax.set_xlabel("Delta Eta")
    ax.set_ylabel("Delta Phi")
    ax.set_title("Angular distance between Higgs")
    outName = "/t3home/gcelotto/ggHH4b/plots/genHiggsEtaPhi.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in", outName)

# NEW PLOT
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    binsEta = np.linspace(-5, 5, 40)
    eta1, eta2 = zip(*etas)
    hist, xedges, yedges = np.histogram2d(eta1, eta2, bins=[binsEta, binsEta])
    im = ax.imshow(hist.T, cmap=plt.cm.viridis, norm=LogNorm(), origin='lower', aspect='auto', extent=(binsEta.min(), binsEta.max(), binsEta.min(), binsEta.max()))
    cbar = fig.colorbar(im, ax=ax, label='Events', shrink=1) 
    ax.set_xlabel("Eta 1")
    ax.set_ylabel("Eta 2")
    outName = "/t3home/gcelotto/ggHH4b/plots/genHiggsEta1Eta2.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)

# newplot
    fig, ax = plt.subplots(1, 2, figsize=(9, 6), constrained_layout=True)
    ax[0].hist(eta1, binsEta, color='blue', label='Eta 1', alpha=0.4)
    ax[0].hist(eta2, binsEta, color='red' , label='Eta 2', alpha=0.4)
    ax[0].set_xlabel("Higgs eta")
    
    ax[1].hist(pt1, binsPt, color='blue', label='Pt 1', alpha=0.4)
    ax[1].hist(pt2, binsPt, color='red' , label='Pt 2', alpha=0.4)
    ax[1].set_xlabel("Higgs pt [GeV]")
    ax[0].legend()
    ax[1].legend()
    
    outName = "/t3home/gcelotto/ggHH4b/plots/hist_higgs_ptEta.png"
    fig.savefig(outName, bbox_inches='tight')
    print("Saving in ", outName)
    return

if __name__ == "__main__":
    main()