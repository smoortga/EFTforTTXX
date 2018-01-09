import os
from argparse import ArgumentParser
import sys
#import matplotlib.pyplot as plt
import numpy as np
from numpy import polyfit, diag, sqrt
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial as P
import numpy as np
import root_numpy as rootnp
from keras.models import load_model
import pickle
import sys
import ROOT
import math
from array import array
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--coupling', default = "cQb1", help='which coupling to use')
parser.add_argument('--TrainingFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/model_checkpoint_save.hdf5", help='path to the directory of SM-like xsec measurement')
parser.add_argument('--ScalerFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/scaler.pkl", help='path to the directory of SM-like xsec measurement')

args = parser.parse_args()

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

def get_Original_nevents(filename,nevents_dict):
    for run,n in nevents_dict.iteritems():
        if run in filename: return n
     
    print "ERROR: couldn't retrieve original number of events from file %s"%filename
    print "Returning -1"
    return -1

def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c  

def GetEventsPassingCut(filelist,variable_name,xmin,xmax,nbins,cut,original_n_events):
    chain = ROOT.TChain("tree")
    for f in filelist:
        chain.Add(args.ValidationDir+"/"+f)
    #print chain.GetEntries()
    tmp_hist = ROOT.TH1D("tmp_hist","",nbins,xmin,xmax)
    chain.Draw(variable_name + " >> tmp_hist")
    #print tmp_hist.GetEntries()
    selected_n_events = tmp_hist.Integral(tmp_hist.FindBin(cut),tmp_hist.GetNbinsX()+1)
    result = float(selected_n_events)/float(original_n_events)
    print "fracion selected: %.2f%%"%(result*100)
    return result


def NN_validate(filename,class_number,cuts,original_n_events = 20000):
    X = rootnp.root2array(args.ValidationDir + "/" + filename[0],"tree")
    X = rootnp.rec2array(X)
    for i in range(len(filename)):
        if i == 0: continue
        X_ = rootnp.root2array(args.ValidationDir + "/" + filename[i],"tree")
        X_ = rootnp.rec2array(X_)
        X = np.concatenate((X,X_))
    model = load_model(args.TrainingFile)
    scaler = pickle.load(open(args.ScalerFile,'r'))
    X = scaler.transform(X)
    if class_number==-1:
        coupling_name = filename[0].split("_")[0]
        #print coupling_name
        coupling_class = classes_dict[coupling_name]
        discr_dict = {}
        for class_n in set(i for j,i in classes_dict.iteritems()):
            discr_dict[class_n] = model.predict(X)[:,class_n]
        #discr = np.asarray([j for jdx,j in enumerate(discr_dict[coupling_class])])
        discr = np.asarray([j/(discr_dict[0][jdx]+discr_dict[coupling_class][jdx]) for jdx,j in enumerate(discr_dict[coupling_class])])
        #discr = np.asarray([(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1])])
    else: discr = model.predict(X)[:,class_number]
    nEvents = len(discr)
    print float(len(discr)), sum_original_n_events, 100*float(len(discr))/float(original_n_events),"%"
    result = {}
    for cut in cuts:
        discr_sel = discr[discr >= cut]
        print "selection efficiency NN cut at %.2f: "%cut ,100*float(len(discr_sel))/float(nEvents)
        result[cut] = (float(len(discr_sel))/float(original_n_events))
    return result
    

def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in val: return -float(val.replace(coupling+"_minus","").replace("p","."))
    else: return float(val.replace(coupling+"_","").replace("p","."))
    
def extract_coupling_string(coupling,val):
    """ extract the string value of a coupling in the directory name"""
    return val.replace(coupling+"_","")
    
def func2(x, a, b, c):
        return (a + b*x + c*x*x)   

def func4(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)    

def func(x, a, b, c):
        return a*(1 + b*x + c*x*x) 
    

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and not 'Limits' in d and args.coupling in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

#print process_dict

validation_files = {}
files = [f for f in os.listdir(args.ValidationDir) if ".root" in f and args.coupling in f]
couplings = set([f.split("_")[0] for f in os.listdir(args.ValidationDir) if ".root" in f and args.coupling in f])
coupl_strengths = set([f.split("_")[1] for f in os.listdir(args.ValidationDir) if ".root" in f and args.coupling in f])
for c in couplings:
    validation_files[c]={}
    for s in coupl_strengths:
        validation_files[c][s]=[]
for f in files:
    coupling_name = f.split("_")[0]
    strength_name = f.split("_")[1]
    validation_files[coupling_name][strength_name].append(f)



# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
coupling_names = []


nevents={"run_01":30000,"run_02":0}

M4b_array = array('d',np.arange(0.1,0.8,0.025))

result_dict = {}
for m in M4b_array:
    result_dict[m]={}
    for order in validation_files[args.coupling]:
        result_dict[m][order]=[]

#print result_dict
dummy_orders = result_dict[0.1].keys()


for o in dummy_orders:
        print "Processing ",o
        ##################
        #
        # Get cross section
        #
        ##################
        if os.path.isfile(args.InputDir + "/" + args.coupling + "/" + args.coupling+"_"+ o + "/cross_section.txt"):
            f_ = open(args.InputDir + "/" + args.coupling + "/" + args.coupling+"_"+ o + "/cross_section.txt", 'r')
            lines = f_.readlines()
            xsec = float(lines[-1].split(" ")[2])
            xsec_error = float(lines[-1].split(" ")[3])
        else: continue
        
        ##################
        #
        # NN output
        #
        ##################
        
        filelist = validation_files[args.coupling][extract_coupling_string(args.coupling,o)]
        if len(filelist) == 0: continue
        sum_original_n_events = 0
        for f in filelist:
            sum_original_n_events += get_Original_nevents(f,nevents)
        frac_passing_evts = NN_validate(filelist,class_number=-1,cuts=M4b_array,original_n_events = sum_original_n_events)
        #result_dict[coupling][extract_coupling_string(coupling,o)] = [xsec,xsec_error,frac_passing_evts, sum_original_n_events]
        #print frac_passing_evts
        for m4b in M4b_array:
            result_dict[m4b][o] = [xsec,xsec_error,frac_passing_evts[m4b], sum_original_n_events]
       
#print result_dict

sm_xsec_frac_error = 0.1
sm_xsec_frac_error_300fb = 0.1

fracerr_measured_sm = sm_xsec_frac_error


lower_limits = array('d')
upper_limits = array('d')
xaxis = array('d')
sm_cross_section_array = array('d')
sm_cross_section_error_array = array('d')

for idx, m4b in enumerate(M4b_array):
    order = result_dict[m4b]
#for idx,m4b,order in result_dict.iteritems():
    coupling_strengths = []
    xsec = []
    xsec_error = []
    sm_values = []
    #nevents = []
    for o,result in order.iteritems():
        if o == "0p0": 
            sm_xsec = result[0]*result[2]
            sm_xsec_error = math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[3]*result[2])/result[3])**2 )
        coupling_strengths.append(extract_coupling(args.coupling,args.coupling+"_"+o))
        xsec.append(result[0]*result[2])
        xsec_error.append( math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[3]*result[2])/result[3])**2 ) )
        #nevents.append(float(lines[-1].split(" ")[4][:-1])) #need to remove \n character in the end
    
   #  SM values
#     sm_xsec = np.mean(np.asarray(sm_values))
#     sm_xsec_error = np.std(np.asarray(sm_values))

    coeff, covmat = curve_fit(func,coupling_strengths,xsec,sigma=xsec_error)#absolute_sigma
    errors = sqrt(diag(covmat))  
    xmin = -20
    xmax = +20
    nsteps = 200
    x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
    y = func(x,*(coeff))
    y_up = func(x,*(coeff + [i for i in errors]))
    y_down = func(x,*(coeff - [i for i in errors])) 
    y_1 = [(i-sm_xsec)**2/((sm_xsec_frac_error*sm_xsec)**2) for i in y_down]
    p_1 = P.fit(x, y_1, 4)
    roots_1 = sorted([i.real for i in (p_1 - 3.84).roots() if i.imag == 0])
    print m4b, roots_1
    while len(roots_1)>2: roots_1 = roots_1[1:-1]
    lower_limits.append(roots_1[0])
    upper_limits.append(roots_1[1])
    xaxis.append(m4b)
    sm_cross_section_array.append(sm_xsec)
    sm_cross_section_error_array.append(sm_xsec_error)
    

# print lower_limits
# print upper_limits
# print xaxis

gr_limit = ROOT.TGraphAsymmErrors( len(xaxis), xaxis, array("d",[0]*len(xaxis)), array("d",[0]*len(xaxis)), array("d",[0]*len(xaxis)),array("d",[abs(i) for i in lower_limits]),array("d",[abs(i) for i in upper_limits]) )
gr_limit.SetFillStyle(3002)
gr_upper_limit = ROOT.TGraph( len(xaxis), xaxis, upper_limits )
gr_lower_limit = ROOT.TGraph( len(xaxis), xaxis, lower_limits )
# gr_upper_limit_300fb = ROOT.TGraph( len(xaxis), xaxis, upper_limits_300fb )
# gr_upper_limit_300fb.SetLineStyle(2)
# gr_lower_limit_300fb = ROOT.TGraph( len(xaxis), xaxis, lower_limits_300fb )
# gr_lower_limit_300fb.SetLineStyle(2)
 
mg = ROOT.TMultiGraph("mg",";NN output;%s [TeV^{-2}]"%convertToLatex(args.coupling))
mg.Add(gr_upper_limit,"l")
mg.Add(gr_lower_limit,"l")
mg.Add(gr_limit,"3")

c = ROOT.TCanvas("c","c",800,800)
uppad = ROOT.TPad("u","u",0.,0.4,1.,1.)
downpad = ROOT.TPad("d","d",0.,0.,1.,0.4)
uppad.Draw()
downpad.Draw()
uppad.cd()
ROOT.gPad.SetMargin(0.15,0.05,0.01,0.1)
mg.Draw("AC")
if "1" in args.coupling:
    mg.SetMinimum(-9)
    mg.SetMaximum(9)
elif "8" in args.coupling:
    mg.SetMinimum(-12)
    mg.SetMaximum(12)
mg.GetYaxis().CenterTitle()
mg.GetYaxis().SetTitleSize(0.06)
mg.GetYaxis().SetTitleOffset(1.05)
mg.GetYaxis().SetLabelSize(0.05)
mg.GetXaxis().SetTitleSize(0.06)
mg.GetXaxis().SetTitleOffset(1.1)
mg.GetXaxis().SetLabelSize(0.05)

l = ROOT.TLegend(0.5,0.8,0.88,0.88)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.SetTextFont(42)
l.SetTextSize(0.055)
l.AddEntry(gr_limit,"95% CL limit @ 300 fb^{-1}","f")
# l.AddEntry(gr_upper_limit_300fb,"95% CL @ 300 fb^{-1}","l")
# l.AddEntry(gr_upper_validity,"EFT validity","l")

l.Draw("same")


#print sm_cross_section_array,sm_cross_section_error_array
downpad.cd()
ROOT.gPad.SetMargin(0.15,0.05,0.25,0.01)
gr_xsec = ROOT.TGraphErrors(len(xaxis), xaxis, array("d",[i*1000 for i in sm_cross_section_array]), array("d",[0]*len(xaxis)),array("d",[i*1000 for i in sm_cross_section_error_array]) )
gr_xsec.SetFillColorAlpha(1,0.7)
mg2 = ROOT.TMultiGraph("mg2",";NN output;cross section [fb]")
mg2.Add(gr_xsec,"3")
mg2.Draw("AC")
mg2.GetYaxis().SetTitleSize(0.095)
mg2.GetYaxis().SetTitleOffset(0.65)
mg2.GetYaxis().SetLabelSize(0.075)
mg2.GetXaxis().CenterTitle()
mg2.GetXaxis().SetTitleSize(0.085)
mg2.GetXaxis().SetTitleOffset(1.3)
mg2.GetXaxis().SetLabelSize(0.085)

l2 = ROOT.TLegend(0.5,0.8,0.88,0.88)
l2.SetFillStyle(0)
l2.SetBorderSize(0)
l2.SetTextFont(42)
l2.SetTextSize(0.09)
l2.AddEntry(gr_xsec,"#splitline{SM cross section after}{NN output > NN^{sel}}","f")
l2.Draw("same")


if not os.path.isdir(args.InputDir+ "/LimitsVsNNOutput"): os.mkdir(args.InputDir+ "/LimitsVsNNOutput")
    
c.SaveAs('%s/LimitsVsNNOutput_%s.png'%(args.InputDir+ "/LimitsVsNNOutput",args.coupling))
c.SaveAs('%s/LimitsVsNNOutput_%s.pdf'%(args.InputDir+ "/LimitsVsNNOutput",args.coupling))   
c.SaveAs('%s/LimitsVsNNOutput_%s.C'%(args.InputDir+ "/LimitsVsNNOutput",args.coupling))   
      
 

