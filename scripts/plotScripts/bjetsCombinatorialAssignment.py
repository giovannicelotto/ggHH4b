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
    fileNames = glob.glob(path+"/ggHH4b*.root")[:10]
    print("%d files to be used" %len(fileNames))
    correct, matchedEvents, totalEntries, correctPair = 0, 0, 0, 0
    higgsMasses = []
    wrongDistances = []
    rightDistances = []
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
            Jet_bReg2018                = branches["Jet_bReg2018"][ev]
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
            p0, p1, p2, p3 = 0.3034,    0.1752,  0.1849,     0.6697 
            scores = p0*(Jet_pt[:jetsToCheck] - np.min(Jet_pt[:jetsToCheck])) / (np.max(Jet_pt[:jetsToCheck]) - np.min(Jet_pt[:jetsToCheck]))
            scores = scores + p1*(1- (abs(Jet_eta[:jetsToCheck]) - np.min(abs(Jet_eta[:jetsToCheck]))) / (np.max(abs(Jet_eta[:jetsToCheck])) - np.min(abs(Jet_eta[:jetsToCheck]))))
            scores = scores + p2*((Jet_qgl[:jetsToCheck] - np.min(Jet_qgl[:jetsToCheck])) / (np.max(Jet_qgl[:jetsToCheck]) - np.min(Jet_qgl[:jetsToCheck])))
            scores = scores + p3*((Jet_btagDeepFlavB[:jetsToCheck] - np.min(Jet_btagDeepFlavB[:jetsToCheck])) / (np.max(Jet_btagDeepFlavB[:jetsToCheck]) - np.min(Jet_btagDeepFlavB[:jetsToCheck])))
            scores = scores + IsTrigJet
            taken = np.argsort(scores)[::-1][:4]


            
            if np.array_equal(np.sort(taken), np.arange(nJet)[m]):
                correct = correct + 1

                j1 = ROOT.TLorentzVector(0.,0.,0.,0.)
                j2 = ROOT.TLorentzVector(0.,0.,0.,0.)
                j3 = ROOT.TLorentzVector(0.,0.,0.,0.)
                j4 = ROOT.TLorentzVector(0.,0.,0.,0.)

                j1.SetPtEtaPhiM(Jet_pt[taken[0]]*Jet_bReg2018[taken[0]], Jet_eta[taken[0]], Jet_phi[taken[0]], Jet_mass[taken[0]])
                j2.SetPtEtaPhiM(Jet_pt[taken[1]]*Jet_bReg2018[taken[1]], Jet_eta[taken[1]], Jet_phi[taken[1]], Jet_mass[taken[1]])
                j3.SetPtEtaPhiM(Jet_pt[taken[2]]*Jet_bReg2018[taken[2]], Jet_eta[taken[2]], Jet_phi[taken[2]], Jet_mass[taken[2]])
                j4.SetPtEtaPhiM(Jet_pt[taken[3]]*Jet_bReg2018[taken[3]], Jet_eta[taken[3]], Jet_phi[taken[3]], Jet_mass[taken[3]])
                def distanceHiggsPlane(m1, m2):
                    x0 = (m1 + m2)/2
                    y0 = x0
                    dist = np.sqrt((m1-x0)**2+(m2-y0)**2)
                    return dist
                
                #print("Invariant masses:")
                #print("1-2 paring mass, mass, distance")
                dist1 = distanceHiggsPlane(m1=(j1+j2).M(), m2=(j3+j4).M())
                #print((j1+j2).M(), (j3+j4).M(), dist1)
                # 1-3 2-4
                dist2 = distanceHiggsPlane(m1=(j1+j3).M(), m2=(j2+j4).M())
                #print("1-3 paring mass, mass, distance")
                #print((j1+j3).M(), (j2+j4).M(), dist2)
                # 1-4 2-3
                dist3 = distanceHiggsPlane(m1=(j1+j4).M(), m2=(j2+j3).M())
                #print("1-4 paring mass, mass, distance")
                #print((j1+j4).M(), (j2+j3).M(), dist3)

                # which is the right assignment:
                for i in taken:
                    for j in taken:
                        if i == j:
                            continue
                        elif GenJet_partonMotherIdx[Jet_genJetIdx[i]] == GenJet_partonMotherIdx[Jet_genJetIdx[j]]:
                            couple1 = [i, j]
                            continue
                couple2 = []
                for i in taken:
                    if i not in couple1:
                        couple2.append(i)
                        
                m1, m2 = np.array([(j1+j2).M(), (j1+j3).M(), (j1+j4).M()]), np.array([(j3+j4).M(), (j2+j4).M(), (j2+j3).M()])
                
                minDistance = np.min((dist1, dist2, dist3))
                minIndex    = np.argmin((dist1, dist2, dist3))
                #newMinIndex = np.argmin((m1+m2)/2)
                #print((m1[minIndex]+m2[minIndex])/2, (m1[newMinIndex]+m2[newMinIndex])/2)

                if abs(minDistance - sorted([dist1, dist2, dist3])[1]) < 25:
                    #print("Here")
                    minIndex    = np.argmin((m1+m2)/2)
                    minDistance = [dist1, dist2, dist3][np.argmin((m1+m2)/2)]

                if minDistance==dist1:
                    if (sorted(couple1)==sorted(list(taken[[0, 1]]))) | (sorted(couple1)==sorted(list(taken[[2, 3]]))):
                        correctPair = correctPair+1
                        rightDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                    else:
                        # wrong
                        wrongDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                        pass
                elif minDistance==dist2:
                    if (sorted(couple1)==sorted(list(taken[[0,2]]))) | (sorted(couple1)==sorted(list(taken[[1,3]]))):
                        correctPair = correctPair+1
                        rightDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                    else:
                        # wrong
                        wrongDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                        pass

                elif minDistance==dist3:
                    if (sorted(couple1)==sorted(list(taken[[0,3]]))) | (sorted(couple1)==sorted(list(taken[[1,2]]))):
                        correctPair = correctPair+1
                        rightDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                    else:
                        # wrong
                        wrongDistances.append(abs(minDistance - sorted([dist1, dist2, dist3])[1]))
                
                    
                higgsMasses.append([m1[minIndex], m2[minIndex]])

    #fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    #higgsMasses = np.array(higgsMasses)
    #bins = np.linspace(0, 300, 100)
    #ax[0].hist(np.clip(higgsMasses[:,0], bins[0], bins[-1]), color='blue', bins= bins)
    #ax[1].hist(np.clip(higgsMasses[:,1], bins[0], bins[-1]), color='blue', bins= bins)
    #fig.savefig("/t3home/gcelotto/ggHH4b/plots/h1mass_h2mass.png", bbox_inches='tight')

    print(matchedEvents/totalEntries)
    print(correct/matchedEvents)
    print(correctPair/correct)

    print("Right distances:", np.mean(rightDistances), np.std(rightDistances))
    print("Wrong distances:", np.mean(wrongDistances), np.std(wrongDistances))
    return

if __name__ == "__main__":
    main()