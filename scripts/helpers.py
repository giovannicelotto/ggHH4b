import glob
import sys
import pandas as pd
import dask.dataframe as dd 
def loadParquet(signalPath, realDataPath, nSignalFiles=-1, nRealDataFiles=1, columns=None):
    signalFileNames = glob.glob(signalPath+"/ggHH*.parquet")
    realDataFileNames = glob.glob(realDataPath+"/Data*.parquet")
    signalFileNames = signalFileNames[:nSignalFiles] if nSignalFiles!=-1 else signalFileNames
    realDataFileNames = realDataFileNames[:nRealDataFiles] if nRealDataFiles!=-1 else realDataFileNames
    print("%d files for MC ggHH4b" %len(signalFileNames))
    print("%d files for realDataFileNames" %len(realDataFileNames))
    
    signal = pd.read_parquet(signalFileNames, columns=columns)
    realData = pd.read_parquet(realDataFileNames, columns=columns)
    return signal, realData

def loadDask(signalPath, realDataPath, nSignalFiles, nRealDataFiles):
    signalFileNames = glob.glob(signalPath+"/*.parquet")
    realDataFileNames = glob.glob(realDataPath+"/*.parquet")
    signalFileNames = signalFileNames[:nSignalFiles] if nSignalFiles!=-1 else signalFileNames
    realDataFileNames = realDataFileNames[:nRealDataFiles] if nRealDataFiles!=-1 else realDataFileNames    

    print("%d files for MC ggHbb" %len(signalFileNames))
    print("%d files for realDataFileNames" %len(realDataFileNames))

    
    try:    
        signal = dd.read_parquet(signalFileNames, blocksize='64 MiB')
        realData = dd.read_parquet(realDataFileNames,  blocksize='64 MiB')
        return signal, realData
    except:
        print("Some of the files might be corrupted. Here is the list:\n")
        for fileName in signalFileNames:
            try:
                df=pd.read_parquet(fileName)
            except:
                print(fileName)
        for fileName in realDataFileNames:
            try:
                df=pd.read_parquet(fileName)
            except:
                print(fileName)
        sys.exit("Exiting the program due to a corrupted files.")



def getXSectionBR():
    xSectionGGHH = 10.517 # fb ## x section at nlo pp->gg->HH->4b  an2019_250_v6 page 7
    return xSectionGGHH