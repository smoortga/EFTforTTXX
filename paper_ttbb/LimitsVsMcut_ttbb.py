import os
from argparse import ArgumentParser
import sys
import ROOT
#import matplotlib.pyplot as plt
import numpy as np
from numpy import polyfit, diag, sqrt
from numpy.polynomial import Polynomial as P
from scipy.optimize import curve_fit
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')
import pickle
from array import array

ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--run', default = "run_02", help='which run to use')
parser.add_argument('--coupling', default = "cQb1", help='which coupling to use')

args = parser.parse_args()



def GetEffMcut(f,mcut):
    return 1.0


def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in val: return -float(val.replace(coupling+"_minus","").replace("p","."))
    else: return float(val.replace(coupling+"_","").replace("p","."))
    
def func(x, a, b, c):
        return a*(1 + b*x + c*x*x)    
 
def convertToLatex(c):
    if len(c) == 4: return c[0].upper()+"_{"+c[1:3]+"}^{"+c[3:4]+"}"
    else: 
        print "couldn't convert"
        return c   

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and args.coupling in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

print process_dict





coupling_names = []

sm_xsec = 0.085 # pb from TOP-16-010
sm_xsec_error = 0.041# pb from TOP-16-010
sm_xsec_frac_error = 0.041/0.085
sm_xsec_frac_error_300fb = 0.1

Mcut_arr = range(600,5101,250)

for coupling,value in process_dict.iteritems():
    if not os.path.isdir(args.InputDir + "/LimitsVsMcut"): os.mkdir(args.InputDir + "/LimitsVsMcut")
    
    coupling_strengths = []
    xsec = []
    xsec_error = []
    nevents = []
    selection_eff_arr = []
    for v in value:
        coupling_strengths.append(extract_coupling(coupling,v))
        f_ = open(args.InputDir + "/" + coupling + "/" + v + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec.append(float(lines[-1].split(" ")[2]))
        xsec_error.append(float(lines[-1].split(" ")[3]))
        nevents.append(float(lines[-1].split(" ")[4][:-1])) #need to remove \n character in the end
        selection_eff = []
        #print args.InputDir + "/" + coupling + "/" + v + "/Events/"+args.run+"/McutDict.pkl"
        delphes_eff_dict = pickle.load(open(args.InputDir + "/" + coupling + "/" + v + "/Events/"+args.run+"/McutDict.pkl","rb"))
        for mcut in Mcut_arr:
            selection_eff.append(delphes_eff_dict[mcut])
        selection_eff_arr.append(selection_eff)
    
    # print xsec
#     print selection_eff_arr
    
    
    
    #**************************
    #
    #   GET LIMITS FOR EACH Mcut
    #
    #**************************
    lower_limits = array('d')
    upper_limits = array('d')
    lower_limits_300fb = array('d')
    upper_limits_300fb = array('d')
    
    for idx,mcut in enumerate(Mcut_arr): 
        eff_tmp = [i[idx] for i in selection_eff_arr]
        xsec_tmp = [i*j for i,j in zip(xsec,eff_tmp)]
        xsec_error_tmp = [i*j for i,j in zip(xsec_error,eff_tmp)]
        #print xsec_tmp
        coeff, covmat = curve_fit(func,coupling_strengths,xsec_tmp,sigma=xsec_error_tmp)#absolute_sigma
        errors = sqrt(diag(covmat))
    
        #print coeff, sqrt(diag(covmat))
        coupling_names.append(coupling)
    
        sm_xsec_tmp = sm_xsec*selection_eff_arr[coupling_strengths.index(0.0)][idx]

        xmin = -40
        xmax = +40
        nsteps = 50
        x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
        #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
        y = func(x,*(coeff))
        y_up = func(x,*(coeff + [i for i in errors]))
        y_down = func(x,*(coeff - [i for i in errors]))    
        y = [(i-sm_xsec_tmp)**2/((sm_xsec_frac_error*sm_xsec_tmp)**2) for i in y]
        ymin = -2
        ymax = 10
        p = P.fit(x, y, 4)
        roots = sorted([i.real for i in (p - 3.84).roots() if i.imag == 0])
        print roots
        while len(roots)>2: roots = roots[1:-1]
        lower_limits.append(roots[0])
        upper_limits.append(roots[1])
        
        
        #******************************************************
        # do the same for optimistic case at 300fb-1
        #******************************************************
        
        y = func(x,*(coeff))
        y = [(i-min(y))**2/((sm_xsec_frac_error_300fb*min(y))**2) for i in y]
        ymin = -2
        ymax = 10
        p = P.fit(x, y, 4)
        roots = sorted([i.real for i in (p - 3.84).roots() if i.imag == 0])
        print roots
        while len(roots)>2: roots = roots[1:-1]
        lower_limits_300fb.append(roots[0])
        upper_limits_300fb.append(roots[1])
        
    
    xaxis = array('d',[i/1000. for i in Mcut_arr])
    
    #print xaxis, upper_limits, lower_limits
    gr_upper_limit = ROOT.TGraph( len(xaxis), xaxis, upper_limits )
    gr_lower_limit = ROOT.TGraph( len(xaxis), xaxis, lower_limits )
    gr_upper_limit_300fb = ROOT.TGraph( len(xaxis), xaxis, upper_limits_300fb )
    gr_upper_limit_300fb.SetLineStyle(2)
    gr_lower_limit_300fb = ROOT.TGraph( len(xaxis), xaxis, lower_limits_300fb )
    gr_lower_limit_300fb.SetLineStyle(2)
    
    # validity limits
    x_validity = array('d',[i/1000. for i in range(min(Mcut_arr),max(Mcut_arr),100)])
    x_validity.append(max(Mcut_arr)/1000.)
    y_upper_validity = array('d',[(16*3.14159*3.14159)/(x*x) for x in x_validity])
    y_lower_validity = array('d',[-(16*3.14159*3.14159)/(x*x) for x in x_validity])
    #print x_validity,y_upper_validity
    
    gr_upper_validity = ROOT.TGraph( len(x_validity), x_validity, y_upper_validity )
    gr_upper_validity.SetLineColor(2)
    #gr_upper_validity.SetLineStyle(2)
    gr_lower_validity = ROOT.TGraph( len(x_validity), x_validity, y_lower_validity )
    gr_lower_validity.SetLineColor(2)
    #gr_Lower_validity.SetLineStyle(2)
    
    mg = ROOT.TMultiGraph("mg",";M_{cut} [TeV];%s [TeV^{-2}]"%convertToLatex(coupling))
    mg.Add(gr_upper_limit)
    mg.Add(gr_lower_limit)
    mg.Add(gr_upper_limit_300fb)
    mg.Add(gr_lower_limit_300fb)
    mg.Add(gr_upper_validity)
    mg.Add(gr_lower_validity)
    
    #ROOT.TGaxis.SetMaxDigits(2)
    c = ROOT.TCanvas("c","c",800,500)
    c.SetMargin(0.15,0.1,0.15,0.1)
    mg.Draw("AC")
    mg.SetMinimum(-90)
    mg.SetMaximum(90)
    mg.GetYaxis().CenterTitle()
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetYaxis().SetTitleOffset(1.05)
    mg.GetYaxis().SetLabelSize(0.05)
    mg.GetXaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetTitleOffset(1.1)
    mg.GetXaxis().SetLabelSize(0.05)
    
    l = ROOT.TLegend(0.55,0.7,0.89,0.88)
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    l.AddEntry(gr_upper_limit,"95% CL @ 2.3 fb^{-1}","l")
    l.AddEntry(gr_upper_limit_300fb,"95% CL @ 300 fb^{-1}","l")
    l.AddEntry(gr_upper_validity,"EFT validity","l")
    
    l.Draw("same")
    
    c.SaveAs('%s/LimitsVsMcut_%s.png'%(args.InputDir+ "/LimitsVsMcut",coupling))
    c.SaveAs('%s/LimitsVsMcut_%s.pdf'%(args.InputDir+ "/LimitsVsMcut",coupling))

    
    # f_.savefig('%s/fitted_xsec/fit_xsec_%s.png'%(args.InputDir,coupling))
#     f_.savefig('%s/fitted_xsec/fit_xsec_%s.pdf'%(args.InputDir,coupling))
#     
#     f_.clf()
    print "done, saved in %s/LimitsVsMcut_%s.png"%(args.InputDir+ "/LimitsVsMcut",coupling)