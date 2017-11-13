import os
from argparse import ArgumentParser
import sys
import multiprocessing
import thread
import subprocess
import time


parser = ArgumentParser()

parser.add_argument('--MGDir', default = os.getcwd(), help='path to the main MadGraph directory (form which you can do ./bin/mg5_aMC)')
parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttcc_inclusive_doubleInsertion", help='path to the directory of all restricted processes')
parser.add_argument('--nevents', type=int, default = 20000, help='number of events to simulate for each process')
parser.add_argument('--n_cpu', type=int, default = 1, help='number of cpu cores to use') #multiprocessing.cpu_count()

args = parser.parse_args()

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd)]

# uncomment if you want to pick only one coupling
# coupling_to_keep = "C8Qd"
# process_dict = {}
# subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and coupling_to_keep in d]
# for subdir in subdirs:
#     process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd)]
# print process_dict


#***************************
#
# create a launch_process.txt card and launch.sh to submit the
# event generation to the cluster easily
#
#****************************
    
for coupling,value in process_dict.iteritems():
    for v in value:
        f_ = open(args.InputDir + "/" + coupling + "/" + v + "/launch_process.txt", 'w')
        f_.write("launch "+args.InputDir + "/" + coupling + "/" + v + "/ \n")
        f_.write(" set nevents %i \n"%args.nevents)
        f_.write(" set nb_core %i \n"%args.n_cpu)
        f_.write("launch "+args.InputDir + "/" + coupling + "/" + v + "/ -i \n")
        f_.write(" print_results --path=%s/cross_section.txt --format=short \n"%(args.InputDir + "/" + coupling + "/" + v))
        f_.close()
        
        if not os.path.isdir(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission"): os.mkdir(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission")
        ff_ = open(args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission/launch.sh", 'w')
        ff_.write("%s/bin/mg5_aMC %s"%(args.MGDir,args.InputDir + "/" + coupling + "/" + v + "/launch_process.txt \n"))
        ff_.close()
        
        basedir = args.InputDir + "/" + coupling + "/" + v + "/localgrid_submission"
        print "qsub -q localgrid -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir)
        os.system("qsub -q localgrid -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir))
        #print "qsub -q express -o %s/script.stdout -e %s/script.stderr -lnodes=1:ppn=%i %s/launch.sh"%(basedir,basedir,args.n_cpu,basedir)
        



print "Done! use 'qstat -u $USER' to monitor samples"



        
#***************************
#
# set nevents to args.nevents in run_card.dat
#
#****************************
# for coupling,value in process_dict.iteritems():
#     for v in value:
#         f_ = open(args.InputDir + "/" + coupling + "/" + v + "/Cards/run_card.dat", 'r')
#         lines_ = f_.readlines()
#         f_.close()
#         ff_ = open(args.InputDir + "/" + coupling + "/" + v + "/Cards/run_card.dat", 'w')
#         for l in lines_:
#             if "nevents" in l: ff_.write("  %i = nevents ! Number of unweighted events requested \n"%args.nevents)
#             else: ff_.write(l)
#         
#         ff_.close()
       