from plot import getShap
import tensorflow as tf
#from tensorflow.compat.v1.keras.backend import get_session
#tensorflow.compat.v1.disable_v2_behavior()
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense
from keras.callbacks import EarlyStopping
from keras.optimizers import Adam
from scipy.integrate import simpson
import sys
sys.path.append('/t3home/gcelotto/ggHH4b/scripts')
from helpers import loadDask, loadParquet, getXSectionBR
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split
import pandas as pd
sys.path.append('/t3home/gcelotto/ggHH4b/scripts/plotScripts/fromFlat/')
def doPlotLoss(fit, outName, earlyStop, patience):

    # "Loss"
    plt.close('all')
    plt.figure(2)
    plt.plot(fit.history['loss'])
    plt.plot(fit.history['val_loss'])
    plt.plot(fit.history['accuracy'])
    plt.plot(fit.history['val_accuracy'])
    plt.title('Model Loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    
    # plt.yscale('log')
    #plt.ylim(ymax = max(min(fit.history['loss']), min(fit.history['val_loss']))*1.4, ymin = min(min(fit.history['loss']),min(fit.history['val_loss']))*0.9)
    plt.ylim(ymin=0, ymax=1)
    ymax = min(fit.history['val_loss'])
    ymin = plt.ylim()[0]
    plt.arrow(x=earlyStop.stopped_epoch-patience-1, y=ymax, dx=0, dy=ymin-ymax, length_includes_head=True, head_length=0.033*(ymin-ymax))
    plt.legend(['Train Loss', 'Val Loss', 'Train Accuracy', 'Val Accuracy'], loc='upper right')
    plt.savefig(outName)
    plt.cla()


def getEllipseCoordinate():
    path = "/t3home/gcelotto/ggHH4b/outputs/xmax_ymax_sigmax_sigmay.npy"
    x_max, y_max, sigma_x, sigma_y = np.load(path)
    return x_max, y_max, sigma_x, sigma_y

def is_point_inside_ellipse(x, y, center_x, center_y, a, b):
    distance = ((x - center_x) / a)**2 + ((y - center_y) / b)**2
    return distance <= 1

def HHclassifier(doTrain, nRealDataFiles):
    hp = {
        'epochs'            : 200,
        'patienceES'        : 10,
        'test_split'        : 0.25,
        'validation_split'  : 0.2,
        'learning_rate'     : 1e-5

          }
    if doTrain:
        realDataPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/Data20181A_2023Nov30/ParkingBPH1/crab_data_Run2018A_part1/231130_120505/flatDataForggHH4b"
        signalPath = "/pnfs/psi.ch/cms/trivcat/store/user/gcelotto/bb_ntuples/nanoaod_ggH/ggHH4b2023Dec20/GluGluToHHTo4B_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/crab_GluGluToHHTo4B/231220_160444/flatData"
        signal, realData = loadParquet(signalPath, realDataPath, nSignalFiles=-1, nRealDataFiles=nRealDataFiles, 
                                       columns=['hh_mass', 'h1_pt','h1_mass',
                                                'h2_pt','h2_Rphi', 'h2_mass',
                                                'j1_pt', 'j2_pt', 'j3_pt', 'j4_pt',
                                                'j1_eta', 'j2_eta', 'j3_eta', 'j4_eta',
                                                'j1_btag', 'j2_btag', 'j3_btag', 'j4_btag',
                                                'dEta_h1h2','tau_h1h2',  
                                                'dPhi_j1j2', 'dR_j1j2','dPhi_j3j4', 'dR_j3j4',
                                                'tau_h1j1', 'tau_h2j3', 'tau_h2j4',
                                                'tau_j1j3', 'dPhi_j1j4', 'tau_j1j4',
                                                'tau_j2j3', 'dPhi_j2j4', 'tau_j2j4',
                                                'ht', 'muon_pt'])

        print("%s Features\n\n"%len(signal.columns))
        #signalRegion
        x_max, y_max, sigma_x, sigma_y = getEllipseCoordinate()
        signal      = signal[is_point_inside_ellipse(signal.h1_mass, signal.h2_mass, x_max, y_max, sigma_x, sigma_y)]
        realData    = realData[is_point_inside_ellipse(realData.h1_mass, realData.h2_mass, x_max, y_max, sigma_x, sigma_y)]
        signal, realData      = signal[signal.j1_pt>20], realData[realData.j1_pt>20]
        signal, realData      = signal[signal.j2_pt>20], realData[realData.j2_pt>20]
        signal, realData      = signal[signal.j3_pt>20], realData[realData.j3_pt>20]
        signal, realData      = signal[signal.j4_pt>20], realData[realData.j4_pt>20]

        signal, realData      = signal[(signal.j1_eta<2.5) & (signal.j1_eta>-2.5)], realData[(realData.j1_eta<2.5) & (realData.j1_eta>-2.5)]
        signal, realData      = signal[(signal.j2_eta<2.5) & (signal.j2_eta>-2.5)], realData[(realData.j2_eta<2.5) & (realData.j2_eta>-2.5)]
        signal, realData      = signal[(signal.j3_eta<2.5) & (signal.j3_eta>-2.5)], realData[(realData.j3_eta<2.5) & (realData.j3_eta>-2.5)]
        signal, realData      = signal[(signal.j4_eta<2.5) & (signal.j4_eta>-2.5)], realData[(realData.j4_eta<2.5) & (realData.j4_eta>-2.5)]


        print("%d events for signal in SR\n%d events for background in SR"%(len(signal), len(realData)))


        y_signal = pd.Series([1]*len(signal))
        y_realData = pd.Series([0]*len(realData))
        X = pd.concat([signal, realData], ignore_index=True)
        Y = pd.concat([y_signal, y_realData], ignore_index=True)

        Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=hp['test_split'], random_state=1999, shuffle=True)
        
        
        Xtrain.to_parquet("/t3home/gcelotto/ggHH4b/outputs/df_NN/Xtrain.parquet")
        Xtest.to_parquet("/t3home/gcelotto/ggHH4b/outputs/df_NN/Xtest.parquet")
        Ytrain.to_csv("/t3home/gcelotto/ggHH4b/outputs/df_NN/Ytrain.csv")
        Ytest.to_csv("/t3home/gcelotto/ggHH4b/outputs/df_NN/Ytest.csv")
    
    Xtrain  = pd.read_parquet("/t3home/gcelotto/ggHH4b/outputs/df_NN/Xtrain.parquet")
    Xtest   = pd.read_parquet("/t3home/gcelotto/ggHH4b/outputs/df_NN/Xtest.parquet")
    Ytrain  = pd.read_csv("/t3home/gcelotto/ggHH4b/outputs/df_NN/Ytrain.csv")['0']
    Ytest   = pd.read_csv("/t3home/gcelotto/ggHH4b/outputs/df_NN/Ytest.csv")['0']
       
    print(len(Xtrain.columns), "COLUMNS")
    # Build the model
    if doTrain:
        model = Sequential()
        model.add(tf.keras.layers.Input(shape = len(signal.columns))) 
        model.add(Dense(units=32, activation='relu', kernel_initializer = tf.keras.initializers.glorot_normal( seed=1999)))
        model.add(Dense(units=16, activation='relu', kernel_initializer = tf.keras.initializers.glorot_normal( seed=1999)))
        model.add(Dense(units=8, activation='relu', kernel_initializer = tf.keras.initializers.glorot_normal( seed=1999)))
        model.add(Dense(units=1, activation='sigmoid', kernel_initializer = tf.keras.initializers.glorot_normal( seed=1999)))
        optimizer = Adam(learning_rate = hp['learning_rate'], beta_1=0.9, beta_2=0.999, epsilon=1e-07, name="Adam") #use_ema=False, bema_momentum=0.99, ema_overwrite_frequency=None, 
        model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
        callbacks=[]
        earlyStop = EarlyStopping(monitor = 'val_loss', patience = hp['patienceES'], verbose = 1, restore_best_weights=True)
        callbacks.append(earlyStop)

        fit = model.fit(Xtrain, Ytrain, epochs=hp['epochs'], callbacks=callbacks, validation_split=hp['validation_split'])
        model.save("/t3home/gcelotto/ggHH4b/scripts/NN/model.h5")
        doPlotLoss(fit=fit, outName="/t3home/gcelotto/ggHH4b/plots/NN/loss.png", earlyStop=earlyStop, patience=hp['patienceES'])

    model = load_model("/t3home/gcelotto/ggHH4b/scripts/NN/model.h5")
    y_predict = model.predict(Xtest)
    yTrain_predict = model.predict(Xtrain)
    thresholds = np.linspace(0, 0.99, 100)
    thresholds = np.concatenate((thresholds, np.linspace(0.99, 1., 1000)))
    
    signal_predictions = y_predict[Ytest==1]
    realData_predictions = y_predict[Ytest==0]
    signalTrain_predictions = yTrain_predict[Ytrain==1]
    print("train predictions", len(signalTrain_predictions))
    print("test predictions", len(signal_predictions))
    realDataTrain_predictions = yTrain_predict[Ytrain==0]
    

    tpr, tpr_train = [], []
    fnr, fnr_train = [], []
    for t in thresholds:
        tpr.append(np.sum(signal_predictions > t)/len(signal_predictions))
        fnr.append(np.sum(realData_predictions > t)/len(realData_predictions))
        tpr_train.append(np.sum(signalTrain_predictions > t)/len(signalTrain_predictions))
        fnr_train.append(np.sum(realDataTrain_predictions > t)/len(realDataTrain_predictions))
    tpr, fnr =np.array(tpr), np.array(fnr)
    tpr_train, fnr_train =np.array(tpr_train), np.array(fnr_train)
    # auc
    auc = -simpson(tpr, fnr)
    auc_train = -simpson(tpr_train, fnr_train)
    print("AUC : ", auc)

    
    fig, ax = plt.subplots(1, 1)#, sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    #fig.subplots_adjust(hspace=0.1)
    ax.plot(fnr, tpr, marker='o', markersize=1, label='test')
    ax.plot(fnr_train, tpr_train, marker='o', markersize=1, label='train')
    ax.plot(thresholds, thresholds, linestyle='dotted', color='green')
    cuts = [0.00001, 0.01, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999]
    bkgRetained = []
    signalRetained = []
    for c in cuts:
        print("BKG retained statistics : ", np.sum(y_predict[Ytest==0]>c))
        bkgRetained.append(np.sum(y_predict[Ytest==0]>c)/len(y_predict[Ytest==0]))
        signalRetained.append(np.sum(y_predict[Ytest==1]>c)/len(y_predict[Ytest==1]))
        print(c, bkgRetained[-1], signalRetained[-1], signalRetained[-1]/np.sqrt(bkgRetained[-1]))
    ax.grid(True)
    ax.set_ylabel("TPR = Signal retained")
    ax.set_xlabel("FNR = Background retained")
    ax.set_xlim(1e-5,1)
    ax.set_ylim(1e-1,1)
    ax.text(x=0.95, y=0.32, s="AUC Test : %.3f"%auc, ha='right')
    ax.text(x=0.95, y=0.28, s="AUC Train: %.3f"%auc_train, ha='right')
    ax.legend()
    #ax.set_yscale('log')
    
    m=fnr>0
    #print(np.max(tpr[m]/np.sqrt(fnr[m])), tpr[np.argmax(tpr[m]/np.sqrt(fnr[m]))], fnr[np.argmax(tpr[m]/np.sqrt(fnr[m]))])
    #    retained.append(tpr[np.argmin(np.array(fnr)[np.array(fnr)>=c])])
    #ax[0].text(x=0.95, y=0.26, s="%s : %.1f%%"%(cuts[0], retained[0]*100), ha='right')
    #ax[0].text(x=0.95, y=0.22, s="%s : %.1f%%"%(cuts[1], retained[1]*100), ha='right')
    #ax[0].text(x=0.95, y=0.18, s="%s : %.1f%%"%(cuts[2], retained[2]*100), ha='right')
    #ax[0].text(x=0.95, y=0.14, s="%s : %.1f%%"%(cuts[3], retained[3]*100), ha='right')
    #ax[0].vlines(ymin=np.zeros(len(retained)), ymax=retained, x=cuts, linestyle='dotted', color='red')
    #ax[0].hlines(xmin=np.zeros(len(retained)), xmax=cuts, y=retained, linestyle='dotted', color='red')

    #m=(fnr>0)
    #ax[1].plot(fnr[m], tpr[m]/np.sqrt(fnr[m]))
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax[1].set_xlabel("FNR")
    #ax[1].set_ylabel("Significance gain")
    fig.savefig("/t3home/gcelotto/ggHH4b/plots/NN/nn_roc.png", bbox_inches='tight')

    fig, ax = plt.subplots(1, 1)
    bins=np.linspace(0, 1, 20)
    sig_test_counts = np.histogram(signal_predictions, bins=bins)[0]
    bkg_test_counts = np.histogram(realData_predictions, bins=bins)[0]
    sig_train_counts = np.histogram(signalTrain_predictions, bins=bins)[0]
    bkg_train_counts = np.histogram(realDataTrain_predictions, bins=bins)[0]
    sig_test_counts, bkg_test_counts, sig_train_counts, bkg_train_counts = sig_test_counts/np.sum(sig_test_counts), bkg_test_counts/np.sum(bkg_test_counts), sig_train_counts/np.sum(sig_train_counts), bkg_train_counts/np.sum(bkg_train_counts)
    ax.hist(bins[:-1], bins=bins, weights=sig_test_counts, alpha=0.3, label='signal test')
    ax.hist(bins[:-1], bins=bins, weights=bkg_test_counts, alpha=.3, label='bkg test')
    ax.scatter((bins[1:]+bins[:-1])/2, sig_train_counts, label='signal train')
    ax.scatter((bins[1:]+bins[:-1])/2, bkg_train_counts, label='bkg train')
    ax.legend(loc='upper center')
    ax.set_yscale('log')
    ax.set_xlabel("NN output")
    fig.savefig("/t3home/gcelotto/ggHH4b/plots/NN/nn_outputs.png", bbox_inches='tight')



    #getShap(Xtest[:1000], model)
    


    return

if __name__ =="__main__":
    doTrain = bool(int(sys.argv[1])) if len(sys.argv)>1 else False
    nRealDataFiles = int(sys.argv[2]) if len(sys.argv)>2 else False
    print("doTrain", doTrain)
    HHclassifier(doTrain=doTrain, nRealDataFiles=nRealDataFiles)