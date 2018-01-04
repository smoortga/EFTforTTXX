import sys
import ROOT
import os
import time
from argparse import ArgumentParser
from array import array
import multiprocessing
import thread
import subprocess
import time
import pickle


ROOT.gROOT.SetBatch(True)

"""
make sure to do the following before running this script!
cdir=$(pwd)
export LD_LIBRARY_PATH=${cdir}:${cdir}/Delphes:$LD_LIBRARY_PATH
"""

def Convert(infilepath, mcut):
    #############
    # input
    #############
    chain = ROOT.TChain("Delphes")
    chain.Add(infilepath)
    
    nEntries = chain.GetEntries()
    
    tmp = ROOT.TH1D("tmp","tmp",100,0,5500)
    chain.Draw("ScalarHT.HT + MissingET.MET >> tmp")
    
    
    counter=tmp.Integral(tmp.GetXaxis().FindBin(0),tmp.GetXaxis().FindBin(mcut))
    
    

    # counter = 0
#     #############
#     # loop
#     #############
#     for evt in range(nEntries):
#         if (evt % int(nEntries/10.) == 0): print"%s: Processing event %i/%i (%.1f %%)"%(infilepath.split("/")[-4] + "_" + infilepath.split("/")[-2] ,evt,nEntries,100*float(evt)/float(nEntries))
#         chain.GetEntry(evt)
# 
#         
#         pass_mcut = True
#         
#         if any([i.HT > mcut for i in chain.ScalarHT]): pass_mcut = False
#         elif any([i.PT> mcut for i in chain.Electron]): pass_mcut = False
#         elif any([i.PT> mcut for i in chain.Muon]): pass_mcut = False
#         elif any([i.PT> mcut for i in chain.Jet]): pass_mcut = False
#         
#         if (pass_mcut): counter += 1
#         
#     del chain

    print "%s: DONE"%(infilepath.split("/")[-4] + "_" + infilepath.split("/")[-2])
    return float(counter)/float(nEntries)

            
      

def main():
    ROOT.gSystem.Load("libDelphes")
    ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
    
    parser = ArgumentParser()
    parser.add_argument('--InputDir', default = "MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
    parser.add_argument('--run', default = "run_02", help='which run to use')
    parser.add_argument('--coupling', default = "cQb1", help='which coupling to use')
    args = parser.parse_args()

    workingdir = os.getcwd()

    #if not os.path.isdir(workingdir+"/CONVERTED_DELPHES_"+args.tag): os.mkdir(workingdir+"/CONVERTED_DELPHES_"+args.tag)
    
    # time_start = time.time()
#     if (args.ncpu < 0 or args.ncpu > multiprocessing.cpu_count()): parallelProcesses = multiprocessing.cpu_count()
#     else: parallelProcesses = args.ncpu
#     p = multiprocessing.Pool(parallelProcesses)
#     print "Using %i parallel processes" %parallelProcesses

    Mcut_arr = range(600,5101,250)
    
    
    coupling_dirs = [i for i in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+i) and not "fitted" in i and args.coupling in i]
    for coupling in coupling_dirs:
        order_dirs = os.listdir(args.InputDir + "/" + coupling)
        for order in order_dirs:
            #if not "Order1" in order and not "Order3" in order: continue
            runs = [i for i in os.listdir(args.InputDir + "/" + coupling + "/" + order + "/Events/") if args.run in i]
            for run in runs:
                root_files =  [i for i in os.listdir(args.InputDir + "/" + coupling + "/" + order + "/Events/"+run) if "delphes_events.root" in i]
                for file in root_files:
                    Mcut_dict={}
                    #print args.InputDir + "/" + coupling + "/" + order + "/Events/"+run+"/"+file
                    input = args.InputDir + "/" + coupling + "/" + order + "/Events/"+run+"/"+file
                    #output = workingdir+"/CONVERTED_DELPHES_"+args.tag+"/"+order+"_"+run+"_"+file
                    #Convert(input,output)
                    for mcut in Mcut_arr:
                        Mcut_dict[mcut] = Convert(input,mcut)
                    print mcut, Mcut_dict
                    pickle.dump(Mcut_dict,open(args.InputDir + "/" + coupling + "/" + order + "/Events/"+args.run+"/McutDict.pkl","wb"))
                    
#     p.close()
#     p.join()
#     
#     print "Total elpased time: %.2f seconds"%(time.time()-time_start)
    
    #Convert("/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed/cQQ1/cQQ1_10p0/Events/run_01/tag_1_delphes_events.root",workingdir+"/CONVERTED_DELPHES_"+args.tag+"/out_test.root")
    
    
    
if __name__ == "__main__":
    main()