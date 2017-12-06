import os
from argparse import ArgumentParser
import sys


#******************************************************
#
# This file takes as an input an example of a param_card.dat file
# and according to a given Black name, it will create restrict_XXX.dat files
# in which all parameters in that block are put to 0 except for a single one
#
#*******************************************************



parser = ArgumentParser()

parser.add_argument('--InputDir', default = os.getcwd(), help='name of the UFO directory')
parser.add_argument('--ExampleParamCard', default = 'param_card.dat', help = 'name of a .dat file in InputDir which contains the block of parameters to restrict')
parser.add_argument('--BlockName', default = 'DIM6', help= 'name of the block that contains the parameters to restrict')
args = parser.parse_args()

parameters_to_scan = [-10.0,-5.0,-1.0,0.0,1.0,5.0,10.0]

#***************************
#
# EXTRACT THE PARAMETERS
#
#****************************

# check if ExampleParamCard exists in InputDir
if args.ExampleParamCard not in os.listdir(args.InputDir):
    print "File %s not found in directory %s"%(args.ExampleParamCard,args.InputDir)
    sys.exit(0)

paramexample = open(args.InputDir+"/"+args.ExampleParamCard, 'r')
line = paramexample.readline()
while not "Block "+args.BlockName in line:
    line = paramexample.readline()
line = paramexample.readline()
params=[]
while len(line.split("#")) == 2:
    coeffiecient = line.split("#")[-1][1:-2] # cut off the " \n" in the end
    params.append(coeffiecient)
    line = paramexample.readline()
    
paramexample.close()

#***************************
#
# FOR EACH PARAMETER: Make a card and put all to zero except the one you need
#
#****************************

summary_restrict_files = open(args.InputDir+"/summary_restrict_files.txt","w")
summary_run_files = open(args.InputDir+"/summary_run_files.txt","w")

for param in params:
    for val in parameters_to_scan:
        parameters_to_ignore = [x for x in params if x != param]
        os.system("cp %s/%s %s/restrict_%s%s.dat"%(args.InputDir,args.ExampleParamCard,args.InputDir,param,str(val).replace(".","p").replace("-","minus")))
        restrictfile = open(args.InputDir+"/restrict_"+param+str(val).replace(".","p").replace("-","minus")+".dat","r")
        lines = restrictfile.readlines()
        restrictfile.close()
    
    
        restrictfile = open(args.InputDir+"/restrict_"+param+str(val).replace(".","p").replace("-","minus")+".dat","w")
        for line in lines:
            parameter_to_ignore_in_line = [p in line for p in parameters_to_ignore]
            if any(parameter_to_ignore_in_line): restrictfile.write(line.replace("1.000000e+00","0.000000e+00"))
            elif param in line: restrictfile.write(line.replace("1.0",str(val)))
            else: restrictfile.write(line)
    
        print "wrote %s"%(args.InputDir+"/restrict_"+param+str(val).replace(".","p").replace("-","minus")+".dat")
        summary_restrict_files.write(args.InputDir+"/restrict_"+param+str(val).replace(".","p").replace("-","minus")+".dat \n")
        
        
        #***************************
        #
        # FOR EACH also write a MACRO for MG5 generation
        #
        #****************************
        
        runfile = open(args.InputDir+"/run_"+param+str(val).replace(".","p").replace("-","minus")+".txt","w")
        if args.InputDir[-1]=="/": name_of_UFO = args.InputDir.split("/")[-2]
        else: name_of_UFO = args.InputDir.split("/")[-1]
        runfile.write("import model ./models/%s-%s \n" %(name_of_UFO,param+str(val).replace(".","p").replace("-","minus")))
        runfile.write("set nb_core 1 \n")
        runfile.write("generate p p > t t~ c c~ QED=0 , (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~); output MODELSCAN_ttcc_DiLepton_ValidationPerCoupling_WithInterference/%s/%s \n"%(param,param+str(val).replace(".","p").replace("-","minus")))
        #runfile.write("generate p p > t t~ c c~ QED=0; output MODELSCAN_ttcc_inclusive_doubleInsertion/%s/%s \n"%(param,param+str(val).replace(".","p").replace("-","minus")))
        runfile.close()
        
        summary_run_files.write(args.InputDir+"/run_"+param+str(val).replace(".","p").replace("-","minus")+".txt \n")
        
        
summary_restrict_files.close()
summary_run_files.close()