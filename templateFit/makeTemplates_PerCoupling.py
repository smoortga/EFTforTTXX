import os
from argparse import ArgumentParser
import sys
import numpy as np
import root_numpy as rootnp
import ROOT
from keras.models import load_model
import pickle
from rootpy.plotting import Hist
#from rootpy.root2array import fill_hist_with_ndarray


def pad(A, length):
    if len(A) > length: return np.asarray(A)
    arr = np.zeros(length)
    arr[:len(A)] = A
    return np.asarray(arr)

def retrieve_discr(input_dir,filename,model,scaler):
    if len(filename) == 0:
        print "ERROR: no files found"
        sys.exit()
    X = rootnp.root2array(input_dir + "/" + filename[0],"tree")
    X = rootnp.rec2array(X)
    for i in range(len(filename)):
        if i == 0: continue
        X_ = rootnp.root2array(input_dir + "/" + filename[i],"tree")
        X_ = rootnp.rec2array(X_)
        X = np.concatenate((X,X_))
    
    X = scaler.transform(X)
    
    discrSM = model.predict(X)[:,0]
    discrLLLL = model.predict(X)[:,1]
    discrRRRR = model.predict(X)[:,2]
    discrLLRR = model.predict(X)[:,3]
    
    #print len(discrSM)
    
    discr_dict = {"SM":np.asarray([]),"LLLL":np.asarray([]),"RRRR":np.asarray([]),"LLRR":np.asarray([])}
    discr_dict["SM"] = discrSM
    discr_dict["LLLL"] = np.asarray([discrLLLL[jdx]/(discrLLLL[jdx]+discrSM[jdx]) for jdx in range(len(discrSM))])
    discr_dict["RRRR"] = np.asarray([discrRRRR[jdx]/(discrRRRR[jdx]+discrSM[jdx]) for jdx in range(len(discrSM))])
    discr_dict["LLRR"] = np.asarray([discrLLRR[jdx]/(discrLLRR[jdx]+discrSM[jdx]) for jdx in range(len(discrSM))])
    
    return discr_dict

def main():

    parser = ArgumentParser()
    parser.add_argument('--TrainingFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/training_output/model_checkpoint_save.hdf5", help='path to training')
    parser.add_argument('--ScalerFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/training_output/scaler.pkl", help='Scaler')
    parser.add_argument('--InputDirSM', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference", help='path to the converted delphes training files')
    parser.add_argument('--InputDirEFT', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ForTemplatesPerCouplings", help='path to the converted delphes training files')
    args = parser.parse_args()
    
    if not os.path.isdir(args.InputDirSM+"/templates/"): os.mkdir(args.InputDirSM+"/templates/")
    f_templates = ROOT.TFile(args.InputDirSM+"/templates/templates.root","RECREATE")
    
    #get training and scaler
    model_ = load_model(args.TrainingFile)
    scaler_ = pickle.load(open(args.ScalerFile,'r'))
    
    # make files_dict
    files_dict = {}
    #SM part
    files_dict["SM"] = []
    fileslist = [f for f in os.listdir(args.InputDirSM) if ".root" in f]
    for f in fileslist:
        if "SM" in f: files_dict["SM"].append(f)
    #EFT part
    files = [f for f in os.listdir(args.InputDirEFT) if ".root" in f]
    couplings = set([f.split("_")[0][:5] if f.split("_")[0][1:3] in ["83","81","13","11"] else f.split("_")[0][:4] for f in os.listdir(args.InputDirEFT) if ".root" in f])
    coupl_strengths = set([f.split("_")[0][5:] if f.split("_")[0][1:3] in ["83","81","13","11"] else f.split("_")[0][4:] for f in os.listdir(args.InputDirEFT) if ".root" in f])
    for c in couplings:
        files_dict[c]={}
        for s in coupl_strengths:
            files_dict[c][s]=[]

    for f in files:
        coupling_name = f.split("_")[0][:5] if f.split("_")[0][1:3] in ["83","81","13","11"] else f.split("_")[0][:4]
        strength_name = f.split("_")[0][5:] if f.split("_")[0][1:3] in ["83","81","13","11"] else f.split("_")[0][4:]
        files_dict[coupling_name][strength_name].append(f)
    
    #print files_dict
    
    # loop over couplings
    for c,arr in files_dict.iteritems():
        print "Processing %s"%c
        if c == "SM":
          discriminators = retrieve_discr(args.InputDirSM,arr,model_,scaler_)
          for output,disc in discriminators.iteritems():
            hist = ROOT.TH1D("h_%s_out%s"%(c,output),";discriminator #frac{P(%s)}{P(%s) + P(SM)};Entries"%(output,output),30,0,1)
            for d in disc:
                hist.Fill(d)
            f_templates.cd()
            hist.Write()
        else:
           for v,files in arr.iteritems():
               print "       %s"%v
               discriminators = retrieve_discr(args.InputDirEFT,files,model_,scaler_)
               for output,disc in discriminators.iteritems():
                  hist = ROOT.TH1D("h_%s_%s_out%s"%(c,v,output),";discriminator #frac{P(%s)}{P(%s) + P(SM)};Entries"%(output,output),30,0,1)
                  for d in disc:
                     hist.Fill(d)
                  f_templates.cd()
                  hist.Write()
                
    
    f_templates.Close()
    

if __name__ == "__main__":
    main()