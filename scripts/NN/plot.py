import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import shap
def getShap(Xtest, model):
        print("Started shap")
        plt.figure()
        max_display = len(Xtest.columns)
        max_display = 32
        explainer = shap.GradientExplainer(model=model, data=Xtest)
        #Compute Shapley values for inX_test[:1000,:]
        print("exapliner started")
        shap_values = explainer.shap_values(np.array(Xtest), nsamples=1000)
         #Generate summary plot
        shap.initjs()
        shap.summary_plot(shap_values, Xtest, plot_type="bar",
                        feature_names=Xtest.columns,
                        max_display=max_display,
                        plot_size=[15.0,0.4*max_display+1.5],
                        class_names=['NN output'],
                        show=False)
        plt.savefig('/t3home/gcelotto/ggHH4b/plots/NN/shap_summary_plot.png')