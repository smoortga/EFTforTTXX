import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--InputDir', default = "FIXME", help='path to the directory of all restricted processes')
args = parser.parse_args()

#args.InputDir = os.getcwd()+"/MODELSCAN_ttcc_inclusive_perOrder_DiLepton"
if not os.path.isdir(args.InputDir):
    print "Could't find directory: %s" %args.InputDir
    sys.exit(1)

coupling_dirs = [i for i in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+i)]
for coupling in coupling_dirs:
    order_dirs = os.listdir(args.InputDir + "/" + coupling)
    for order in order_dirs:
        events_dir = args.InputDir + "/" + coupling + "/" + order + "/Events/"
        runs = [i for i in os.listdir(events_dir) if "run_" in i]
        for run in runs:
            files_to_delete = [i for i in os.listdir(args.InputDir + "/" + coupling + "/" + order + "/Events/"+run) if ".gz" in i]
            for file in files_to_delete:
                path_to_delete = args.InputDir + "/" + coupling + "/" + order + "/Events/"+run+"/"+file
                print "removing file: %s of size %.2f MB"%(path_to_delete,os.path.getsize(path_to_delete)/float(1000000))
                os.system("rm %s"%path_to_delete)
                
