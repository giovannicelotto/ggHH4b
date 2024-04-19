import numpy as np
import uproot
import glob
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
hep.style.use("CMS")
def main(saveAndLoad, nFiles):
    if saveAndLoad:
        path = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000"
        fileNames = glob.glob(path+"/ggHH4b*.root")[:nFiles]
        print("%d files to be used" %len(fileNames))
        matchedEvents, totalEntries = 0, 0
        listOfDf=[]
        for fileName in fileNames:
            f = uproot.open(fileName)
            tree = f['Events']
            branches = tree.arrays()
            maxEntries = tree.num_entries 
            totalEntries = totalEntries + maxEntries
            print("Entries : %d" %maxEntries)

            for ev in  range(maxEntries):


                if (ev%(int(maxEntries/20))==0):
                    sys.stdout.write('\r')
                    # the exact output you're looking for:
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
                df = {
            'h1j1_pt' : [],         'h1j1_eta' : [],         'h1j1_phi' : [],         'h1j1_mass' : [],         'h1j1_btag' : [],         'h1j1_index' : [],         'h1j1_genPt' : [],         'h1j1_genEta' : [],         'h1j1_genPhi' : [],         'h1j1_genMass' : [], 

            'h1j2_pt' : [],         'h1j2_eta' : [],         'h1j2_phi' : [],         'h1j2_mass' : [],         'h1j2_btag' : [],         'h1j2_index' : [],         'h1j2_genPt' : [],         'h1j2_genEta' : [],         'h1j2_genPhi' : [],         'h1j2_genMass' : [], 

            'h2j1_pt' : [],         'h2j1_eta' : [],         'h2j1_phi' : [],         'h2j1_mass' : [],         'h2j1_btag' : [],         'h2j1_index' : [],         'h2j1_genPt' : [],         'h2j1_genEta' : [],         'h2j1_genPhi' : [],         'h2j1_genMass' : [], 

            'h2j2_pt' : [],         'h2j2_eta' : [],         'h2j2_phi' : [],         'h2j2_mass' : [],         'h2j2_btag' : [],         'h2j2_index' : [],         'h2j2_genPt' : [],         'h2j2_genEta' : [],         'h2j2_genPhi' : [],         'h2j2_genMass' : [], 

            'otherJets_pt' : [],         'otherJets_eta' : [],         'otherJets_phi' : [],         'otherJets_mass' : [],         'otherJets_btag' : []    
            }
                matchedEvents = matchedEvents + 1

                h1_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]==GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
                h2_GenJetIdx = [idx for idx in Jet_genJetIdx[m] if GenJet_partonMotherIdx[idx]!=GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
    # np.arange(nJet)[m] are the recojets matched
                hh_RecoJetsIdx = np.arange(nJet)[m]
                h1_RecoJetsIdx = [idx for idx in hh_RecoJetsIdx if GenJet_partonMotherIdx[Jet_genJetIdx[idx]]==GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
                h2_RecoJetsIdx = [idx for idx in hh_RecoJetsIdx if GenJet_partonMotherIdx[Jet_genJetIdx[idx]]!=GenJet_partonMotherIdx[Jet_genJetIdx[m]][0]]
                df['h1j1_pt'].append(np.float32(Jet_pt[h1_RecoJetsIdx[0]]))
                df['h1j1_eta'].append(np.float32(Jet_eta[h1_RecoJetsIdx[0]]))
                df['h1j1_phi'].append(np.float32(Jet_phi[h1_RecoJetsIdx[0]]))
                df['h1j1_mass'].append(np.float32(Jet_mass[h1_RecoJetsIdx[0]]))
                df['h1j1_btag'].append(np.float32(Jet_btagDeepFlavB[h1_RecoJetsIdx[0]]))
                df['h1j1_index'].append(np.float32(h1_RecoJetsIdx[0]))
                df['h1j1_genPt'].append(np.float32(GenJet_pt[h1_GenJetIdx[0]]))
                df['h1j1_genEta'].append(np.float32(GenJet_eta[h1_GenJetIdx[0]]))
                df['h1j1_genPhi'].append(np.float32(GenJet_phi[h1_GenJetIdx[0]]))
                df['h1j1_genMass'].append(np.float32(GenJet_mass[h1_GenJetIdx[0]]))

                df['h1j2_pt'].append(np.float32(Jet_pt[h1_RecoJetsIdx[1]]))
                df['h1j2_eta'].append(np.float32(Jet_eta[h1_RecoJetsIdx[1]]))
                df['h1j2_phi'].append(np.float32(Jet_phi[h1_RecoJetsIdx[1]]))
                df['h1j2_mass'].append(np.float32(Jet_mass[h1_RecoJetsIdx[1]]))
                df['h1j2_btag'].append(np.float32(Jet_btagDeepFlavB[h1_RecoJetsIdx[1]]))
                df['h1j2_index'].append(np.float32(h1_RecoJetsIdx[1]))
                df['h1j2_genPt'].append(np.float32(GenJet_pt[h1_GenJetIdx[1]]))
                df['h1j2_genEta'].append(np.float32(GenJet_eta[h1_GenJetIdx[1]]))
                df['h1j2_genPhi'].append(np.float32(GenJet_phi[h1_GenJetIdx[1]]))
                df['h1j2_genMass'].append(np.float32(GenJet_mass[h1_GenJetIdx[1]]))

                df['h2j1_pt'].append(np.float32(Jet_pt[h2_RecoJetsIdx[0]]))
                df['h2j1_eta'].append(np.float32(Jet_eta[h2_RecoJetsIdx[0]]))
                df['h2j1_phi'].append(np.float32(Jet_phi[h2_RecoJetsIdx[0]]))
                df['h2j1_mass'].append(np.float32(Jet_mass[h2_RecoJetsIdx[0]]))
                df['h2j1_btag'].append(np.float32(Jet_btagDeepFlavB[h2_RecoJetsIdx[0]]))
                df['h2j1_index'].append(np.float32(h2_RecoJetsIdx[0]))
                df['h2j1_genPt'].append(np.float32(GenJet_pt[h2_GenJetIdx[0]]))
                df['h2j1_genEta'].append(np.float32(GenJet_eta[h2_GenJetIdx[0]]))
                df['h2j1_genPhi'].append(np.float32(GenJet_phi[h2_GenJetIdx[0]]))
                df['h2j1_genMass'].append(np.float32(GenJet_mass[h2_GenJetIdx[0]]))

                df['h2j2_pt'].append(np.float32(Jet_pt[h2_RecoJetsIdx[1]]))
                df['h2j2_eta'].append(np.float32(Jet_eta[h2_RecoJetsIdx[1]]))
                df['h2j2_phi'].append(np.float32(Jet_phi[h2_RecoJetsIdx[1]]))
                df['h2j2_mass'].append(np.float32(Jet_mass[h2_RecoJetsIdx[1]]))
                df['h2j2_btag'].append(np.float32(Jet_btagDeepFlavB[h2_RecoJetsIdx[1]]))
                df['h2j2_index'].append(np.float32(h2_RecoJetsIdx[1]))
                df['h2j2_genPt'].append(np.float32(GenJet_pt[h2_GenJetIdx[1]]))
                df['h2j2_genEta'].append(np.float32(GenJet_eta[h2_GenJetIdx[1]]))
                df['h2j2_genPhi'].append(np.float32(GenJet_phi[h2_GenJetIdx[1]]))
                df['h2j2_genMass'].append(np.float32(GenJet_mass[h2_GenJetIdx[1]]))

                bkgJetsIdx = np.arange(nJet)[~ m]
                for idx in bkgJetsIdx:
                    df['otherJets_pt'].append(np.float32(Jet_pt[idx]))
                    df['otherJets_eta'].append(np.float32(Jet_eta[idx]))
                    df['otherJets_phi'].append(np.float32(Jet_phi[idx]))
                    df['otherJets_mass'].append(np.float32(Jet_mass[idx]))
                    df['otherJets_btag'].append(np.float32(Jet_btagDeepFlavB[idx]))
                listOfDf.append(df)
        df=pd.DataFrame(listOfDf)
        df.to_parquet("signalAndOtherJets.parquet")
    df = pd.read_parquet("signalAndOtherJets.parquet")
    
    #print(np.concatenate([df.h1j1_pt, df.h1j2_pt, df.h2j1_pt, df.h2j2_pt]))
    
    fig, ax = plt.subplots(2, 2, constrained_layout=True)
    bins_pt = np.linspace(0, 200, 30)
    bins_eta = np.linspace(-5, 5, 30)
    bins_mass = np.linspace(0, 20, 30)
    bins_btag = np.linspace(0, 1, 10)
    counts_pt = np.histogram(np.concatenate([df.h1j1_pt, df.h1j2_pt, df.h2j1_pt, df.h2j2_pt]), bins=bins_pt)[0]
    counts_eta = np.histogram(np.concatenate([df.h1j1_eta, df.h1j2_eta, df.h2j1_eta, df.h2j2_eta]), bins=bins_eta)[0]
    counts_mass = np.histogram(np.concatenate([df.h1j1_mass, df.h1j2_mass, df.h2j1_mass, df.h2j2_mass]), bins=bins_mass)[0]
    counts_btag = np.histogram(np.concatenate([df.h1j1_btag, df.h1j2_btag, df.h2j1_btag, df.h2j2_btag]), bins=bins_btag)[0]
    counts_pt   = counts_pt/np.sum(counts_pt)
    counts_eta  = counts_eta/np.sum(counts_eta)
    counts_mass = counts_mass/np.sum(counts_mass)
    counts_btag = counts_btag/np.sum(counts_btag)
    ax[0,0].hist(bins_pt[:-1], bins=bins_pt, weights=counts_pt, histtype=u'step', color='blue', label='Signal')
    ax[0,1].hist(bins_eta[:-1], bins=bins_eta, weights=counts_eta, histtype=u'step', color='blue', label='Signal')
    ax[1,0].hist(bins_mass[:-1], bins=bins_mass, weights=counts_mass, histtype=u'step', color='blue', label='Signal')
    ax[1,1].hist(bins_btag[:-1], bins=bins_btag, weights=counts_btag, histtype=u'step', color='blue', label='Signal')
    
    
    counts_pt   = np.histogram(np.clip(df.otherJets_pt.explode('otherJets_pt'), bins_pt[0], bins_pt[-1]), bins=bins_pt)[0]
    counts_eta  = np.histogram(np.clip(df.otherJets_eta.explode('otherJets_eta'), bins_eta[0], bins_eta[-1]), bins=bins_eta)[0]
    counts_mass = np.histogram(np.clip(df.otherJets_mass.explode('otherJets_mass'), bins_mass[0], bins_mass[-1]), bins=bins_mass)[0]
    counts_btag = np.histogram(np.clip(df.otherJets_btag.explode('otherJets_btag'), bins_btag[0], bins_btag[-1]), bins=bins_btag)[0]
    counts_pt   = counts_pt/np.sum(counts_pt)
    counts_eta  = counts_eta/np.sum(counts_eta)
    counts_mass = counts_mass/np.sum(counts_mass)
    counts_btag = counts_btag/np.sum(counts_btag)
    ax[0, 0].hist(bins_pt[:-1], bins=bins_pt, weights=counts_pt, histtype=u'step', color='red', label='Other Jets')
    ax[0, 1].hist(bins_eta[:-1], bins=bins_eta, weights=counts_eta, histtype=u'step', color='red', label='Other Jets')
    ax[1, 0].hist(bins_mass[:-1], bins=bins_mass, weights=counts_mass, histtype=u'step', color='red', label='Other Jets')
    ax[1, 1].hist(bins_btag[:-1], bins=bins_btag, weights=counts_btag, histtype=u'step', color='red', label='Other Jets')

    ax[0, 0].set_xlabel("Jet pt [GeV]")
    ax[0, 1].set_xlabel("Jet eta ")
    ax[1, 0].set_xlabel("Jet mass [GeV]")
    ax[1, 1].set_xlabel("Jet btag")

    ax[0, 0].legend()
    ax[0, 1].legend()
    ax[1, 0].legend()
    ax[1, 1].legend()
    fig.savefig("SignalVsOtherJets.png", bbox_inches='tight')


# determine njets to look at
    df['max_index'] = df[['h1j1_index', 'h1j2_index', 'h2j1_index', 'h2j2_index']].max(axis=1)

    countsVsMax = []
    for i in range(3, 15):
        
        m = df['max_index'] < i
        countsVsMax.append(np.sum(m))
    countsVsMax=np.array(countsVsMax)
    print(countsVsMax)
    countsVsMax=countsVsMax/len(df)
    print(countsVsMax)
    fig, ax = plt.subplots(1, 1)
    bins=np.arange(3, 15)
    print(bins, countsVsMax)
    ax.hist(bins, bins=bins+0.5, weights=countsVsMax, histtype=u'step', color='blue')
    for i in range(len(bins)):
        ax.text(x=bins[i]-0.5, y=countsVsMax[i], s="%.1f"%(countsVsMax[i]*100), fontsize=20)
    fig.savefig("/t3home/gcelotto/ggHH4b/maxIdx.png")
    
    sys.exit("Exiting")



if __name__ == "__main__":
    if len(sys.argv)==2:
        saveAndLoad=bool(int(sys.argv[1]))
        main(saveAndLoad, nFiles=None)
    if len(sys.argv)>2:
        saveAndLoad=bool(int(sys.argv[1]))
        nFiles=int(sys.argv[2])
        main(saveAndLoad, nFiles)