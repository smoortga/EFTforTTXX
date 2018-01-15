import os
from argparse import ArgumentParser
import sys
import multiprocessing
import thread
import subprocess
import time
import sys




def extract_coupling(coupling,val,handed="l"):
    """ extract the numeric value of a coupling in the directory name"""
    if handed == "l": 
        value = val.split("_")[0] + "_" + val.split("_")[1]
        coupl = coupling.split("_")[0]
    elif handed == "r":
        value = val.split("_")[2] + "_" + val.split("_")[3]
        coupl = coupling.split("_")[1]
    else: 
        print "ERROR"
        return -9999999
    if "minus" in value: return -float(value.replace(coupl+"_minus","").replace("p",".").replace("_EFT","").replace("_SM",""))
    else: return float(value.replace(coupl+"_","").replace("p",".").replace("_EFT","").replace("_SM",""))


parser = ArgumentParser()

parser.add_argument('--MGDir', default = os.getcwd(), help='path to the main MadGraph directory (form which you can do ./bin/mg5_aMC)')
parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_Marginalized", help='path to the directory of all restricted processes')
parser.add_argument('--nevents', type=int, default = 30000, help='number of events to simulate for each process')
parser.add_argument('--n_cpu', type=int, default = 1, help='number of cpu cores to use') #multiprocessing.cpu_count()

args = parser.parse_args()


fourH_names = [
    "cQQ1",
    "cQQ8",
    "cQt1",
    "cQb1",
    #"ctt1", # does not create diagrams at LO EFT for ttbb
    "ctb1",
    "cQt8",
    "cQb8",
    "ctb8"
]


process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd)]



#***************************
#
# create a launch_process.txt card and launch.sh to submit the
# event generation to the cluster easily
#
#****************************
    
for coupling,value in process_dict.iteritems():
    for v in value:
        lefthanded = coupling.split("_")[0]
        righthanded = coupling.split("_")[1]
        left_handed_value = extract_coupling(coupling,v,"l")
        right_handed_value = extract_coupling(coupling,v,"r")
        #if not "Order1" in v and not "Order3" in v: continue
        workingdir = os.getcwd()
        os.chdir(args.InputDir + "/" + coupling + "/" + v + "/Source/")
        os.system("make clean > /dev/null")
        os.chdir(workingdir)
    
        f_ = open(args.InputDir + "/" + coupling + "/" + v + "/launch_process.txt", 'w')
        f_.write("set nb_core %i \n"%args.n_cpu)
        f_.write("launch "+args.InputDir + "/" + coupling + "/" + v + "/ \n")
        f_.write(" shower=PYTHIA6 \n")
        f_.write(" detector=DELPHES \n")
        f_.write(" done \n")
        f_.write(" /user/smoortga/Analysis/MG5_aMC_v2_6_0/Delphes/cards/delphes_card_EFT.dat \n")
        f_.write(" set nevents %i \n"%args.nevents)
        for c_ in fourH_names:
            if c_ != lefthanded and c_ != righthanded: f_.write(" set %s 0.0 \n"%c_)
            elif c_ == lefthanded: f_.write(" set %s %s \n"%(c_,str(left_handed_value)))
            elif c_ == righthanded: f_.write(" set %s %s \n"%(c_,str(right_handed_value)))
        f_.write(" set ctt1 0.0 \n")
        f_.write(" set ptl 20 \n")
        f_.write(" set ptb 20 \n")
        f_.write(" set etaj 2.5 \n")
        f_.write(" set etab 2.5 \n")
        f_.write(" set drjl 0.0 \n")
        f_.write(" set drbl 0.0 \n")
        f_.write(" done \n")
        f_.write("launch "+args.InputDir + "/" + coupling + "/" + v + "/ -i \n")
        f_.write(" print_results --path=%s/cross_section.txt --format=short \n"%(args.InputDir + "/" + coupling + "/" + v))
        f_.close()
        
        if not os.path.isdir(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission"): os.mkdir(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission")
        ff_ = open(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission/launch.sh", 'w')
        ff_.write("#!/bin/bash \n")
        ff_.write("pwd=$PWD \n")  
        ff_.write("source $VO_CMS_SW_DIR/cmsset_default.sh \n")                                                                                                                                                           
        ff_.write("cd /storage_mnt/storage/user/smoortga/CMSSW_8_0_21/src \n")                                                                                                                                                          
        ff_.write("eval `scram runtime -sh` \n")                                                                                                                                           
        ff_.write("cd $pwd \n") 
        ff_.write("%s/bin/mg5_aMC %s"%(args.MGDir,args.InputDir + "/" + coupling + "/" + v + "/launch_process.txt \n"))
        ff_.close()
        
        basedir = args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission"
        print "qsub -q localgrid -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir)
        os.system("qsub -q localgrid -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir))
        #print "qsub -q express -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir)
        



print "Done! use 'qstat -u $USER' to monitor samples"
print "use 'for j in $(qselect -u $USER);do timeout 3 qdel -a $j;done' to delete all your jobs"


