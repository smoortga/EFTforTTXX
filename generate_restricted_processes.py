import os
from argparse import ArgumentParser
import sys
import multiprocessing
import thread
import subprocess
import time

#******************************************************
#
# This file takes as an input a file like "summary_run_files.txt" which is output of
# the write_restriction_cards.py script. It runs MG on these predefined commands, generating the process
# that is defined in each run_XXXX.txt and outputting this to a specific directory
#
#*******************************************************




parser = ArgumentParser()

parser.add_argument('--InputDir', default = os.getcwd(), help='path to the directory of input text file containing the run_XXX.txt paths')
parser.add_argument('--InputFile', default = "summary_run_files.txt", help='name of input text file containing the run_XXX.txt paths')
parser.add_argument('--n_cpu', type=int, default = 4, help='number of cpu cores to use') #multiprocessing.cpu_count()

args = parser.parse_args()





# small helper functions for assync processing


def runMG(textfile):
    os.system("./bin/mg5_aMC %s"%textfile)
    #write_launch_file(textfile)




# check if InputFile exists in InputDir
if args.InputFile not in os.listdir(args.InputDir):
    print "File %s not found in directory %s"%(args.InputFile,args.InputDir)
    sys.exit(0)

infile = open(args.InputDir+"/"+args.InputFile, 'r')
lines = [l.split(" ")[0] for l in infile.readlines()] # cut off the " \n" in the end

# run in parallel n_cpu generations in MG    
time_start = time.time()
parallelProcesses = args.n_cpu
p = multiprocessing.Pool(parallelProcesses)
print "Using %i parallel processes" %parallelProcesses
print "%i jobs to run" %(len(lines))

for run_card_path in lines:
    if not os.path.isfile(run_card_path): 
        print "%s not found, skipping..."%run_card_path
        continue
    p.apply_async(runMG, args = (run_card_path,))

p.close()
p.join()

print "Total elpased time: %.2f seconds"%(time.time()-time_start)

