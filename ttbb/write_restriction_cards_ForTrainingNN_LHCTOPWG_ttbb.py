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

#orders_to_scan = [0,1,2,3,4]
param_train_value = 3.0
param_train_string = "3p0"

#***************************
#
# EXTRACT THE PARAMETERS
#
#****************************

# check if ExampleParamCard exists in InputDir
if args.ExampleParamCard not in os.listdir(args.InputDir):
    print "File %s not found in directory %s"%(args.ExampleParamCard,args.InputDir)
    sys.exit(0)

twoHtwoL_names = [
    "cQq83",
    "cQq81",
    "cQu8",
    "cQd8",
    "ctq8",
    "ctu8",
    "ctd8",
    "cQq13",
    "cQq11",
    "cQu1",
    "cQd1",
    "ctq1",
    "ctu1",
    "ctd1"
]

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

paramexample = open(args.InputDir+"/"+args.ExampleParamCard, 'r')
line = paramexample.readline()
while not "Block "+args.BlockName in line:
    line = paramexample.readline()
line = paramexample.readline()
params=[]
params_to_delete=[]
while len(line.split("#")) == 2:
    coefficient = line.split("#")[-1][1:-2] # cut off the " \n" in the end
    if coefficient not in fourH_names: 
        line = paramexample.readline()
        params_to_delete.append(coefficient)
        continue
    params.append(coefficient)
    line = paramexample.readline()
    
paramexample.close()

print params
#sys.exit(1)

#***************************
#
# FOR EACH PARAMETER: Make a card and put all to zero except the one you need
#
#****************************

summary_restrict_files = open(args.InputDir+"/summary_restrict_files.txt","w")
summary_run_files = open(args.InputDir+"/summary_run_files.txt","w")

for param in params:

    parameters_to_ignore = [x for x in params if x != param]+params_to_delete
    os.system("cp %s/%s %s/restrict_%s_3p0.dat"%(args.InputDir,args.ExampleParamCard,args.InputDir,param))
    restrictfile = open(args.InputDir+"restrict_"+param+"_3p0"+".dat","r")
    lines = restrictfile.readlines()
    restrictfile.close()


    restrictfile = open(args.InputDir+"/restrict_"+param+"_3p0"+".dat","w")
    for line in lines:
        if "Lambda" in line: 
            restrictfile.write(line)
            continue
        parameter_to_ignore_in_line = [p in line for p in parameters_to_ignore]
        #if any(parameter_to_ignore_in_line): restrictfile.write(line.replace("1.000000e+00","0.000000e+00"))
        if any(parameter_to_ignore_in_line): restrictfile.write(line.replace("1.000000e+00","0.000000e+00"))
        elif param in line: restrictfile.write(line.replace("1.0",str(param_train_value)))
        else: restrictfile.write(line)
        #else: restrictfile.write(line.replace("1.000000e+00","3.000000e+00"))

    print "wrote %s"%(args.InputDir+"/restrict_"+param+"_3p0"+".dat")
    summary_restrict_files.write(args.InputDir+"/restrict_"+param+"_3p0"+".dat \n")
    
    
    #***************************
    #
    # write pure EFT  MACRO for MG5 generation
    #
    #****************************
    
    runfile = open(args.InputDir+"/run_"+param+"_3p0_EFT.txt","w")
    if args.InputDir[-1]=="/": name_of_UFO = args.InputDir.split("/")[-2]
    else: name_of_UFO = args.InputDir.split("/")[-1]
    runfile.write("import model ./models/%s-%s \n" %(name_of_UFO,param+"_3p0"))
    runfile.write("set nb_core 1 \n")
    runfile.write("define p = p b b~ \n")
    runfile.write("define j = p \n")
    #runfile.write("generate p p > t t~ c c~ , (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~) QED=0; output MODELSCAN_ttcc_dilepton/%s/%s \n"%(param,param+str(val).replace(".","p").replace("-","minus")))
    runfile.write("generate p p > t t~ b b~ QED=0 DIM6=1 QCD=2, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~); output MODELSCAN_ttbb_DiLepton_Training_WithSpinCorrelations/%s/%s \n"%(param,param+"_3p0_EFT"))
    runfile.close()
    
    print "wrote %s"%(args.InputDir+"/run_"+param+"_3p0_EFT.txt")
    summary_run_files.write(args.InputDir+"/run_"+param+"_3p0_EFT.txt \n")
    
    #***************************
    #
    # write SM MACRO for MG5 generation
    #
    #****************************
    
    runfile = open(args.InputDir+"/run_"+param+"_3p0_SM.txt","w")
    if args.InputDir[-1]=="/": name_of_UFO = args.InputDir.split("/")[-2]
    else: name_of_UFO = args.InputDir.split("/")[-1]
    runfile.write("import model ./models/%s-%s \n" %(name_of_UFO,param+"_3p0"))
    runfile.write("set nb_core 1 \n")
    # for 5-flavour scheme
    runfile.write("define p = p b b~ \n")
    runfile.write("define j = p \n")
    #runfile.write("generate p p > t t~ c c~ , (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~) QED=0; output MODELSCAN_ttcc_dilepton/%s/%s \n"%(param,param+str(val).replace(".","p").replace("-","minus")))
    runfile.write("generate p p > t t~ b b~ QED=0 DIM6=0 QCD=4, (t > w+ b, w+ > l+ vl), (t~ > w- b~, w- > l- vl~); output MODELSCAN_ttbb_DiLepton_Training_WithSpinCorrelations/%s/%s \n"%(param,param+"_3p0_SM"))
    runfile.close()
    
    print "wrote %s"%(args.InputDir+"/run_"+param+"_3p0_SM.txt")
    summary_run_files.write(args.InputDir+"/run_"+param+"_3p0_SM.txt \n")
        
        
summary_restrict_files.close()
summary_run_files.close()