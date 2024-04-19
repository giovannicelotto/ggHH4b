import uproot
import glob
import numpy as np
import pandas as pd
import ROOT
'''
    take all the MINIAOD files of ggHH4b and get the number of entries, save them.
    take all the NANOAOD files of ggHH4b and get the number of entries
    take all the flat    files of ggHH4b and get the number of entries
'''

def getMiniEntries():
    fileNames = [
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/AFC17129-F999-5D4D-B5D1-A16C3549EE53.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/79557878-F2C1-AD4B-8F69-60A7F12EF423.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/58E08C9A-11E3-3345-ADCE-67976E6C6F13.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/3D1E7373-2E13-C546-AE25-E069779FC2AC.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/27BF2AE5-45A8-A84A-8565-0B10E7FCC4C0.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/110BD57F-3B7E-CD47-97DC-7EB29036299F.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/4C229D15-4048-E54A-B08F-3C6E9AC21F07.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/D30ED9B5-7EB7-2945-939D-A78619B955B5.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/31A4FB1C-4E74-9A47-8161-A535097E930D.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/E324F3E9-1F5B-F448-A22A-A5C8A3C5D184.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/1A4BE7EE-6DEF-0642-80E0-2367BAA60239.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/07816ED8-708F-1C42-B9EA-C87C39394C3F.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/A0C18898-DDFD-334F-9B3E-2F6F2CC99746.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/F5DE7AB6-187B-D44A-B6A9-895430B66575.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/3D66F53B-88E1-0245-9BE4-940D222F05EA.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/EF7C22DF-D418-1242-86D2-6F2EDF3FEB6C.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/A0557F3F-0E82-8C41-9C87-E8181B81684D.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/D83C4E19-49DD-BD48-96A0-6217A63ABB3A.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/024CC802-B97B-E246-88A6-86EAFB1BC356.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/41C57715-E13D-9F41-AA28-FF76D0FCC78C.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/BC708D5D-1598-264D-A3AD-D63ABF80BC3F.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/1D36D9EA-B722-D443-981B-58A62C41D820.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/5D978182-0A61-B341-8DEB-6C76613F3E22.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/EF85CA8F-E7B4-6F4C-BBD8-F9CD44301400.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/15110B3A-6DAD-2448-A1F9-A3D415069ED4.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/EA8DC637-A975-DB4E-8EBA-874C60D4702A.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/E0C9C79E-DAAD-494A-A5F7-8B9DADDCCC7E.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/B9E07207-0EDF-CF4C-B37A-019E73166CAA.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/0AF89BD7-1FB6-2A4D-A669-38FADBE07E71.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/667C2AF5-4E15-F74F-9757-780C81882521.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/39AA1B84-612F-614A-B6A9-712374ED5523.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/C05F0FD1-87C0-E346-96A1-848CCF91E335.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/2B7ECC82-DAC0-9E43-8775-A50F7AF6A659.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/B50DA4F2-2489-FC4E-A3B4-E2CB8FFDCC1E.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/563A092A-C43C-C64B-A86C-C687E1137972.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/D62B569D-017D-D94A-A055-413A805EECB3.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/2DD6A39E-63AD-2D4A-A4F7-DE96936D570F.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/E2856A82-807A-6E4E-8623-3DA7ABE0DA6E.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/21461F62-FD9F-6E4E-ACFD-0435B4949EB2.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/6C61870A-6EBB-FC47-B922-5460BE294170.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/881A5E82-7F47-3E42-A02E-BF415F94C408.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/FE2106C7-935D-7843-85C5-60EAD31CB97E.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/95403BC2-124F-074C-AE6E-7523FBFBB0D4.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/B7641697-E22C-4843-9C19-371B42135AD5.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/330075D9-999F-2148-AF39-71BF08275015.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/048BA098-081F-FF46-8321-4233E838B201.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/F09DCA90-5548-BE40-A034-273691EF33B1.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/11C28A68-853A-AD4E-9204-01130D258332.root',
'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2810000/3993CC83-5358-BB42-938E-BADCB34EE511.root',

]
    totalMiniEntries = 0
    for fileName in fileNames:
        f = ROOT.TFile.Open(fileName)
        tree = f.Get('Events')
        maxEntries = tree.GetEntries()
        totalMiniEntries += maxEntries
        print("%d/%d\n\t\t"%(fileNames.index(fileName)+1, len(fileNames)), totalMiniEntries)
    return totalMiniEntries

def getNanoEntries():
    fileNames = glob.glob("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/0000/*.root")
    print("Number of Flat Data : ", len(fileNames))
    totalNanoEntries = 0
    for fileName in fileNames:
        f = uproot.open(fileName)
        tree = f['Events']
        maxEntries = tree.num_entries 
        totalNanoEntries += maxEntries
        print("%d/%d\n\t\t"%(fileNames.index(fileName)+1, len(fileNames)), totalNanoEntries)
    return totalNanoEntries

def getFlatEntries():
    fileNames = glob.glob("/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData/*.parquet")
    print("Number of Flat Data : ", len(fileNames))
    totalFlatEntries = 0
    for fileName in fileNames:
        f = pd.read_parquet(fileName, columns=['h1_mass'])
        maxEntries = len(f)
        totalFlatEntries += maxEntries
        print("%d/%d\n\t\t"%(fileNames.index(fileName)+1, len(fileNames)), totalFlatEntries)
    return totalFlatEntries



if __name__=="__main__":
    mini = getMiniEntries()
    print("Saving mini into /t3home/gcelotto/ggHH4b/outputs/N_mini.npy ...")
    np.save("/t3home/gcelotto/ggHH4b/outputs/N_mini.npy", mini)

    nano = getNanoEntries()
    print("Saving nano into /t3home/gcelotto/ggHH4b/outputs/N_nano.npy ...")
    np.save("/t3home/gcelotto/ggHH4b/outputs/N_nano.npy", nano)
    
    flat = getFlatEntries()

    efficiencyMC = flat/mini
    
    print("Mini : \t%d"%mini)
    print("Nano : \t%d"%nano)
    print("Flat : \t%d"%flat)
    print("Efficiency_MC  : \t%.10f"%(flat/mini))