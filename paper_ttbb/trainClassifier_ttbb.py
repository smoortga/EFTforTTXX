from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TLegend, gROOT, gStyle
import sys
import ROOT
import os
import time
from argparse import ArgumentParser
from array import array
from math import *
import numpy as np
from collections import Counter
import root_numpy as rootnp
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.utils import np_utils
from keras.callbacks import EarlyStopping, ModelCheckpoint, LearningRateScheduler
from keras.optimizers import SGD,Adam
from keras.regularizers import l1, l2
#from keras.utils.visualize_util import plot
from numpy.lib.recfunctions import stack_arrays
from sklearn.preprocessing import StandardScaler
from keras.models import load_model
from sklearn.metrics import roc_curve,roc_auc_score
from sklearn.cross_validation import train_test_split
import pickle
from rootpy.plotting import Hist

#from rootpy.plotting import Hist2D

from sklearn.neural_network import MLPClassifier



def makeROC(fpr, tpr, thresholds,AUC,outfile,signal_label, background_label):
	
	c = TCanvas("c","c",700,600)
	ROOT.gPad.SetMargin(0.15,0.07,0.15,0.05)
	ROOT.gPad.SetLogy(0)
	ROOT.gPad.SetGrid(1,1)
	ROOT.gStyle.SetGridColor(17)
	
	roc = ROOT.TGraph(len(fpr),tpr,fpr)
	
	roc.SetLineColor(2)
	roc.SetLineWidth(2)
	roc.SetTitle(";Signal efficiency (%s); Background efficiency (%s)"%(signal_label, background_label))
	roc.GetXaxis().SetTitleOffset(1.4)
	roc.GetXaxis().SetTitleSize(0.045)
	roc.GetYaxis().SetTitleOffset(1.4)
	roc.GetYaxis().SetTitleSize(0.045)
	roc.GetXaxis().SetRangeUser(0,1)
	roc.GetYaxis().SetRangeUser(0.000,1)
	roc.Draw("AL")
	
	ROOT.gStyle.SetTextFont(42)
	t = ROOT.TPaveText(0.2,0.84,0.4,0.94,"NBNDC")
	t.SetTextAlign(11)
	t.SetFillStyle(0)
	t.SetBorderSize(0)
	t.AddText('AUC = %.3f'%AUC)
	t.Draw('same')
	
	c.SaveAs(outfile)

def makeDiscr(discr_dict,outfile,xtitle="discriminator"):
    c = ROOT.TCanvas("c","c",800,500)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetMargin(0.15,0.1,0.2,0.1)
    #ROOT.gPad.SetLogy(1)
    #ROOT.gPad.SetGrid(1,1)
    ROOT.gStyle.SetGridColor(17)
    l = TLegend(0.17,0.75,0.88,0.88)
    l.SetTextSize(0.055)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetNColumns(2)
    
    colors = [2,4,8,ROOT.kCyan+2]
    counter = 0
    for leg,discr in discr_dict.iteritems():
        a = Hist(30, 0, 1)
        #fill_hist_with_ndarray(a, discr)
        a.fill_array(discr)
        a.SetLineColor(colors[counter])
        a.SetLineWidth(2)
        a.GetXaxis().SetTitle(xtitle)
        a.GetXaxis().SetLabelSize(0.05)
        a.GetXaxis().SetTitleSize(0.05)
        a.GetXaxis().SetTitleOffset(1.45)
        a.GetYaxis().SetTitle("a.u.")
        a.GetYaxis().SetTickSize(0)
        a.GetYaxis().SetLabelSize(0)
        a.GetYaxis().SetTitleSize(0.06)
        a.GetYaxis().SetTitleOffset(0.9)
        a.Scale(1./a.Integral())
        #a.GetYaxis().SetRangeUser(0.00001,100)
        a.GetYaxis().SetRangeUser(0,0.2)
        if counter == 0: a.draw("hist")
        else: a.draw("same hist")  
        l.AddEntry(a,leg,"l")
        counter += 1
        
    l.Draw("same") 
    
    c.SaveAs(outfile)
    

def drawTrainingCurve(input,output):
    hist = pickle.load(open(input,"rb"))
    tr_acc = hist["acc"]
    tr_loss = hist["loss"]
    val_acc = hist["val_acc"]
    val_loss = hist["val_loss"]
    epochs = range(len(tr_acc))

    plt.figure(1)
    plt.subplot(211)
    plt.plot(epochs, tr_acc,label="training")
    plt.plot(epochs, val_acc, label="validation")
    plt.legend(loc='best')
    plt.grid(True)
    #plt.xlabel("number of epochs")
    plt.ylabel("accuracy")
    plt.subplot(212)
    plt.plot(epochs, tr_loss, label="training")
    plt.plot(epochs, val_loss, label="validation")
    plt.legend(loc='best')
    plt.grid(True)
    plt.xlabel("number of epochs")
    plt.ylabel("loss")
    plt.savefig(output)

  
def make_model(input_dim, nb_classes, nb_hidden_layers = 1, nb_neurons = 50,momentum_sgd = 0.8, init_learning_rate_sgd = 0.0005, dropout =0.1,nb_epoch = 100, batch_size=128):
    #batch_size = 128
    #nb_epoch = args.n_epochs

    #prepare the optimizer 
    decay_sgd = init_learning_rate_sgd/float(5*nb_epoch) if nb_epoch !=0 else 0.0001
    sgd = SGD(lr=init_learning_rate_sgd, decay=decay_sgd, momentum=momentum_sgd, nesterov=True)


    model = Sequential()
    model.add(Dense(nb_neurons ,input_shape= input_dim))
    model.add(Activation('relu'))
    for x in range ( nb_hidden_layers ):
            model.add(Dense(nb_neurons))
            model.add(Activation('relu'))
            model.add(Dropout(dropout))
    # model.add(Dense(nb_neurons))
#     model.add(Activation('relu'))
    #model.add(Dropout(dropout))
    model.add(Dense(nb_classes))
    model.add(Activation('softmax'))

    model.compile(loss='categorical_crossentropy', optimizer=sgd, metrics=['accuracy'])
 
    return model

def main():
    
    gROOT.SetBatch(1)

    parser = ArgumentParser()
    parser.add_argument('--nepoch', type=int, default=100,help='number of epochs to run the training for')
    parser.add_argument('--TrainingFile', default = "", help='path to training')
    parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed", help='path to the converted delphes files')
    #parser.add_argument('--tag', default=time.strftime("%a%d%b%Y_%Hh%Mm%Ss"),help='name of output directory')
    args = parser.parse_args()
    
    
    #couplings = ["SM", "C83Qq", "C81Qq", "C13Qq", "C11Qq", "C8td", "C8tu", "C1td", "C1tu", "C8Qd", "C8Qu", "C8tq", "C1Qd", "C1Qu", "C1tq"]
    #couplings = ["SM","C8td", "C8tu"]
    # classes_dict = { #name:class_number
#         'SM': 0,
#         #LLLL
#         'C11Qq': 1,
#         'C81Qq': 2,
#         'C13Qq': 1,
#         'C83Qq': 2,
#         #RRRR
#         'C1tu': 1,
#         'C8tu': 2,
#         'C1td': 1,
#         'C8td': 2,
#          #LLRR
#         'C8tq': 2,
#         'C8Qu': 2,
#         'C8Qd': 2,
#         'C1tq': 1,
#         'C1Qu': 1,
#         'C1Qd': 1    
#     }
    classes_dict = { #name:class_number
            'SM': 0,
            #LLLL
            "cQQ1": 1,
            "cQQ8": 1,
            #LLRR
            "cQt1": 2,
            "cQb1": 1,
            "cQt8": 2,
            "cQb8": 1,
            #RRRR
            "ctb1": 2,
            "ctb8": 2
        }
#     classes_dict = { #name:class_number
#         'SM': 0,
#         #LLLL
#         #'C11Qq': 1,
#         #'C81Qq': 1,
#         #'C13Qq': 1,
#         #'C83Qq': 1,
#         #RRRR
#         'C1tu': 1,
#         'C8tu': 2,
#         'C1td': 1,
#         'C8td': 2,
#          #LLRR
#         #'C8tq': 1,
#         #'C8Qu': 1,
#         #'C8Qd': 2,
#         #'C1tq': 1,
#         #'C1Qu': 1,
#         #'C1Qd': 2    
#     }
    couplings = classes_dict.keys()
    print couplings
    
    nb_classes = len(set(i for j,i in classes_dict.iteritems()))
    
    files = [i for i in os.listdir(args.InputDir) if ".root" in i]
    
    chain_dict = {}
    for c in couplings:
        chain_dict.update({c:TChain("tree")})
    #SM_chain = TChain("")
    #C1tu_chain = TChain("")


    for f in files:
        if "SM" in f: chain_dict["SM"].Add(args.InputDir + "/" + f)
        elif "EFT" in f:
            coupling_name = f.split("_")[0]#[:-3]
            if coupling_name in couplings: chain_dict[coupling_name].Add(args.InputDir + "/" + f)

    branchnames = [i.GetName() for i in chain_dict["SM"].GetListOfBranches()]
    print branchnames, len(branchnames)

    data_dict = {}
    X_SM_ = rootnp.tree2array(chain_dict["SM"])
    X_SM_ = rootnp.rec2array(X_SM_)
    data_dict.update({"SM":[X_SM_,np.zeros(len(X_SM_))]})
    #counter = 1
    for name, chain in chain_dict.iteritems():
        if name == "SM": continue
        X_ = rootnp.tree2array(chain)
        X_ = rootnp.rec2array(X_)
        y_= np.asarray([classes_dict[name]]*len(X_))
        data_dict.update({name:[X_,y_]})
        #counter += 1
    
    # make sure that all classes have the same number of events
    nb_events_per_class = [0]*nb_classes
    for name, data in data_dict.iteritems():
        nb_events_per_class[classes_dict[name]] += len(data[0])
    smallest_nb_events_per_class = min(nb_events_per_class)
    ratio_per_class = [smallest_nb_events_per_class/float(i) for i in nb_events_per_class]
    
    #shortest_length = min([len(i[0]) for name, i in data_dict.iteritems() ])
    data_dict["SM"][0] = data_dict["SM"][0][0:int(ratio_per_class[classes_dict["SM"]]*len(data_dict["SM"][0]))]
    data_dict["SM"][1] = data_dict["SM"][1][0:int(ratio_per_class[classes_dict["SM"]]*len(data_dict["SM"][1]))]
    X = data_dict["SM"][0]
    y = data_dict["SM"][1]
    for name, data in data_dict.iteritems():
        if name == "SM": continue
        X = np.concatenate((X,data[0][0:int(ratio_per_class[classes_dict[name]]*len(data[0]))]))
        y = np.concatenate((y,data[1][0:int(ratio_per_class[classes_dict[name]]*len(data[1]))]))
    Y = np_utils.to_categorical(y.astype(int), nb_classes)
    
    scaler = StandardScaler()
    scaler.fit(X)
    if not os.path.isdir(args.InputDir + "/training_output"): os.mkdir(args.InputDir + "/training_output")
    pickle.dump(scaler,open(args.InputDir + "/training_output/scaler.pkl",'wb'))
        
    X_train, X_test , y_train, y_test, Y_train, Y_test = train_test_split(X, y, Y, test_size=0.2)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    
    
    
    if args.TrainingFile == "":
        model = make_model(X_train.shape[1:],nb_classes, nb_epoch = args.nepoch)
        print model.summary()
    
        batch_size = 128
        if not os.path.isdir(args.InputDir + "/training_output"): os.mkdir(args.InputDir + "/training_output")
        train_history = model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=args.nepoch, validation_data=(X_test, Y_test), callbacks = [ModelCheckpoint(args.InputDir + "/training_output/model_checkpoint_save.hdf5")], shuffle=True,verbose=1)
                  #sample_weight=sw)

        pickle.dump(train_history.history,open(args.InputDir + "/training_output/loss_and_acc.pkl",'wb'))
    
    else:
        model = load_model(args.TrainingFile)
    
    drawTrainingCurve(args.InputDir + "/training_output/loss_and_acc.pkl",args.InputDir + "/training_output/training_curve.pdf")
    
    discr_dict = {}
    for class_number in set(i for j,i in classes_dict.iteritems()):
        discr_dict[class_number] = model.predict(X_test)[:,class_number]
    
    
    #SM vs EFT
    SM_discr = [(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==0]
    EFT_discr = [(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==1 or y_test[jdx] == 2]
    #tR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==2]
    #LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==3]
    fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(EFT_discr)))),np.concatenate((SM_discr,EFT_discr)))
    AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(EFT_discr)))),np.concatenate((SM_discr,EFT_discr)))
    makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsEFT.pdf","EFT","SM")
    makeDiscr({"EFT":EFT_discr,"SM":SM_discr},args.InputDir + "/training_output/discr_SMvsEFT.pdf","discriminator P(t_{L}) + P(t_{R})")
    
    # SM vs tR
#     SM_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==0]
#     tL_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==1]
#     tR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==2]
#     LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==3]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(tR_discr)))),np.concatenate((SM_discr,tR_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(tR_discr)))),np.concatenate((SM_discr,tR_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvstR.png","t_{R}","SM")
#     makeDiscr({"SM":SM_discr,"t_{L}":tL_discr, "t_{R}":tR_discr},args.InputDir + "/training_output/discr_SMvstR.png","discriminator #frac{P(t_{R})}{P(t_{R}) + P(SM)}")
#     
    #tL vs tR
    #SM_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==0]
    tL_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==1]
    tR_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==2]
    #LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==3]
    fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(tR_discr)),np.ones(len(tL_discr)))),np.concatenate((tR_discr,tL_discr)))
    AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(tR_discr)),np.ones(len(tL_discr)))),np.concatenate((tR_discr,tL_discr)))
    makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_tLvstR.pdf","t_{L}","t_{R}")
    makeDiscr({"#splitline{EFT with}{left-handed top}":tL_discr, "#splitline{EFT with}{right-handed top}":tR_discr},args.InputDir + "/training_output/discr_tLvstR.pdf","discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}")
    

#     SM vs singlet
#     SM_discr = [j for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==0]
#     singlet_discr = [j for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==1]
#     Octet_discr = [j for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==2]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(singlet_discr)))),np.concatenate((SM_discr,singlet_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(singlet_discr)))),np.concatenate((SM_discr,singlet_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsSinglet.png","singlet","SM")
#     makeDiscr({"SM":SM_discr,"color singlet":singlet_discr, "color octet":Octet_discr},args.InputDir + "/training_output/discr_SMvsSinglet.png","discriminator P(singlet)")
#     
#     SM vs Octet
#     SM_discr = [j for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==0]
#     singlet_discr = [j for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==1]
#     Octet_discr = [j for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==2]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(Octet_discr)))),np.concatenate((SM_discr,Octet_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(Octet_discr)))),np.concatenate((SM_discr,Octet_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsOctet.png","octet","SM")
#     makeDiscr({"SM":SM_discr,"color singlet":singlet_discr, "color octet":Octet_discr},args.InputDir + "/training_output/discr_SMvsOctet.png","discriminator P(octet)")
#     
#     Singlet vs Octet
#     SM_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==0]
#     Singlet_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==1]
#     Octet_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==2]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.ones(len(Singlet_discr)),np.zeros(len(Octet_discr)))),np.concatenate((Singlet_discr,Octet_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.ones(len(Singlet_discr)),np.zeros(len(Octet_discr)))),np.concatenate((Singlet_discr,Octet_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SingletvsOctet.png","singlet","octet")
#     makeDiscr({"SM":SM_discr,"color singlet":Singlet_discr, "color octet":Octet_discr},args.InputDir + "/training_output/discr_SingletvsOctet.png","discriminator P(octet)/(P(singlet) + P(octet))")
#     
    #SM vs LLLL
    # SM_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==0]
#     LLLL_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==1]
#     RRRR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==2]
#     LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[1][jdx]) for jdx,j in enumerate(discr_dict[1]) if y_test[jdx] ==3]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(LLLL_discr)))),np.concatenate((SM_discr,LLLL_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(LLLL_discr)))),np.concatenate((SM_discr,LLLL_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsLLLL.png","LLLL","SM")
#     makeDiscr({"SM":SM_discr,"LLLL":LLLL_discr, "RRRR":RRRR_discr,"LLRR":LLRR_discr},args.InputDir + "/training_output/discr_SMvsLLLL.png","discriminator #frac{P(LLLL)}{P(LLLL) + P(SM)}")
#     
#     #SM vs RRRR
#     SM_discr = [j/(discr_dict[0][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==0]
#     LLLL_discr = [j/(discr_dict[0][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==1]
#     RRRR_discr = [j/(discr_dict[0][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==2]
#     LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==3]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(RRRR_discr)))),np.concatenate((SM_discr,RRRR_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(RRRR_discr)))),np.concatenate((SM_discr,RRRR_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsRRRR.png","RRRR","SM")
#     makeDiscr({"SM":SM_discr,"LLLL":LLLL_discr, "RRRR":RRRR_discr,"LLRR":LLRR_discr},args.InputDir + "/training_output/discr_SMvsRRRR.png","discriminator #frac{P(RRRR)}{P(RRRR) + P(SM)}")
#     
#     #SM vs LLRR
#     SM_discr = [j/(discr_dict[0][jdx]+discr_dict[3][jdx]) for jdx,j in enumerate(discr_dict[3]) if y_test[jdx] ==0]
#     LLLL_discr = [j/(discr_dict[0][jdx]+discr_dict[3][jdx]) for jdx,j in enumerate(discr_dict[3]) if y_test[jdx] ==1]
#     RRRR_discr = [j/(discr_dict[0][jdx]+discr_dict[3][jdx]) for jdx,j in enumerate(discr_dict[3]) if y_test[jdx] ==2]
#     LLRR_discr = [j/(discr_dict[0][jdx]+discr_dict[3][jdx]) for jdx,j in enumerate(discr_dict[3]) if y_test[jdx] ==3]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(LLRR_discr)))),np.concatenate((SM_discr,LLRR_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(SM_discr)),np.ones(len(LLRR_discr)))),np.concatenate((SM_discr,LLRR_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_SMvsLLRR.png","LLRR","SM")
#     makeDiscr({"SM":SM_discr,"LLLL":LLLL_discr, "RRRR":RRRR_discr,"LLRR":LLRR_discr},args.InputDir + "/training_output/discr_SMvsLLRR.png","discriminator #frac{P(LLRR)}{P(LLRR) + P(SM)}")
#     
#     #LLLL vs RRRR
#     SM_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==0]
#     LLLL_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==1]
#     RRRR_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==2]
#     LLRR_discr = [j/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[2]) if y_test[jdx] ==3]
#     fpr, tpr, thres = roc_curve(np.concatenate((np.zeros(len(LLLL_discr)),np.ones(len(RRRR_discr)))),np.concatenate((LLLL_discr,RRRR_discr)))
#     AUC = 1-roc_auc_score(np.concatenate((np.zeros(len(LLLL_discr)),np.ones(len(RRRR_discr)))),np.concatenate((LLLL_discr,RRRR_discr)))
#     makeROC(fpr, tpr, thres,AUC,args.InputDir + "/training_output/roc_LLLLvsRRRR.png","RRRR","LLLL")
#     makeDiscr({"SM":SM_discr,"LLLL":LLLL_discr, "RRRR":RRRR_discr,"LLRR":LLRR_discr},args.InputDir + "/training_output/discr_LLLLvsRRRR.png","discriminator #frac{P(LLLL)}{P(LLLL) + P(RRRR)}")
#     





    # X_sig = rootnp.tree2array(C1tu_chain)
#     X_sig = rootnp.rec2array(X_sig)
#     X_bkg = rootnp.tree2array(SM_chain)
#     X_bkg = rootnp.rec2array(X_bkg)
#     X = np.concatenate((X_sig,X_bkg))
#     y = np.concatenate((np.ones(len(X_sig)),np.zeros(len(X_bkg))))
#     X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
# 
#     #mlp_parameters = {'early_stopping':list([True,False]),'activation':list(['tanh','relu']), 'hidden_layer_sizes':list([(5,10),(10,15),(20,50)]), 'algorithm':list(['adam']), 'alpha':list([0.0001,0.00005]), 'tol':list([0.00001]), 'learning_rate_init':list([0.001,0.005,0.0005])}
#     #mlp_parameters = {'activation':list(['relu']), 'hidden_layer_sizes':list([(50,50,50)]), 'algorithm':list(['adam']), 'alpha':list([0.00001]), 'tol':list([0.00001]), 'learning_rate_init':list([0.001])}
#     mlp_clf = MLPClassifier(max_iter = 100, activation='relu', hidden_layer_sizes=(5,10), alpha=0.001,learning_rate_init=0.001, learning_rate = 'adaptive',verbose=1,tol=0.00001) #learning_rate = 'adaptive'
#     mlp_clf.fit(X_train,y_train)
# 
#     mlp_disc = mlp_clf.predict_proba(X_test)[:,1]
#     mlp_fpr, mlp_tpr, mlp_thresholds = roc_curve(y_test, mlp_disc)
#     AUC = 1-roc_auc_score(y_test,mlp_disc)
# 	
#     makeROC(mlp_fpr, mlp_tpr, mlp_thresholds,AUC,"./roc.png")
	
    
if __name__ == "__main__":
    main()