import os
from argparse import ArgumentParser
import sys
#import matplotlib.pyplot as plt
import numpy as np
from numpy import polyfit, diag, sqrt
from numpy.polynomial import Polynomial as P
from scipy.optimize import curve_fit
import ROOT
from array import array
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')

ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
#parser.add_argument('--SMxsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODEL_ttcc_dilepton_SMRestrictedMassAndCabibbo", help='path to the directory of SM-like xsec measurement')

args = parser.parse_args()




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
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and not 'Limits' in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

#print process_dict

# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
coupling_names = []

sm_xsec = 0.085 # pb from TOP-16-010
sm_xsec_error = 0.041# pb from TOP-16-010
sm_xsec_frac_error = 0.041/0.085
sm_xsec_frac_error_300fb = 0.1



# output files for limits
if not os.path.isdir(args.InputDir + "/fitted_xsec"): os.mkdir(args.InputDir + "/fitted_xsec")

limits_file = open('%s/fitted_xsec/limits_ttbb_xsec_CMS.txt'%(args.InputDir),"w")
limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*sm_xsec_frac_error))

limits_file_300fb = open('%s/fitted_xsec/limits_ttbb_xsec_CMS_300fb.txt'%(args.InputDir),"w")
limits_file_300fb.write("fractional error on SM measurement: %i %% \n"%int(100*sm_xsec_frac_error_300fb))

for coupling,value in process_dict.iteritems():
    
    
    coupling_strengths = []
    xsec = []
    xsec_error = []
    nevents = []
    for v in value:
        coupling_strengths.append(extract_coupling(coupling,v))
        f_ = open(args.InputDir + "/" + coupling + "/" + v + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec.append(float(lines[-1].split(" ")[2]))
        xsec_error.append(float(lines[-1].split(" ")[3]))
        nevents.append(float(lines[-1].split(" ")[4][:-1])) #need to remove \n character in the end
    
    coeff, covmat = curve_fit(func,coupling_strengths,xsec,sigma=xsec_error)#absolute_sigma
    errors = sqrt(diag(covmat))
    
    coupling_names.append(coupling)
    p0.append(coeff[0])
    p0e.append(errors[0])
    p1.append(coeff[1])
    p1e.append(errors[1])
    p2.append(coeff[2])
    p2e.append(errors[2])
    
    #**************************
    #
    #   PLOTTING
    #
    #**************************
    
    gr_data = ROOT.TGraphErrors(len(coupling_strengths),array('d',coupling_strengths),array('d',xsec),array('d',np.zeros(len(coupling_strengths))),array('d',xsec_error))
    gr_data.SetLineWidth(1)
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(4)
    xmin = -40
    xmax = +40
    nsteps = 200
    x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
    y = func(x,*(coeff))
    y_up = func(x,*(coeff + [i for i in errors]))
    y_down = func(x,*(coeff - [i for i in errors]))
    gr_fitted_xsec = ROOT.TGraphErrors(len(x),array('d',x),array('d',y),array('d',np.zeros(len(x))),array('d',[abs(i-j) for i,j in zip(y,y_up)]))
    gr_fitted_xsec.SetMarkerSize(0)
    gr_fitted_xsec.SetFillColor(2)
    gr_fitted_xsec.SetLineColor(2)
    gr_fitted_xsec.SetLineWidth(2)
    
    # CMS result @ 2.3 fb-1
    gr_CMS = ROOT.TGraphErrors(len(x),array('d',x),array('d',[sm_xsec]*len(x)),array('d',np.zeros(len(x))),array('d',[2*sm_xsec_frac_error*sm_xsec]*len(x)))
    gr_CMS.SetFillColorAlpha(3,0.5)
    gr_CMS.SetLineWidth(2)
    # CMS result @ 300 fb-1
    gr_CMS_300fb = ROOT.TGraphErrors(len(x),array('d',x),array('d',[min(y)]*len(x)),array('d',np.zeros(len(x))),array('d',[2*sm_xsec_frac_error_300fb*min(y)]*len(x)))
    gr_CMS_300fb.SetFillColorAlpha(5,0.7)
    gr_CMS_300fb.SetLineWidth(2)
    gr_CMS_300fb.SetLineStyle(2)
    
    mg = ROOT.TMultiGraph("mg",";;cross section [pb]")
    mg.Add(gr_CMS,"3")
    mg.Add(gr_CMS_300fb,"3")
    mg.Add(gr_fitted_xsec,"c3")
    mg.Add(gr_data,"pe")
    
    c = ROOT.TCanvas("c","c",800,700)
    uppad = ROOT.TPad("u","u",0.,0.4,1.,1.)
    downpad = ROOT.TPad("d","d",0.,0.,1.,0.4)
    uppad.Draw()
    downpad.Draw()
    uppad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.01,0.1)
    #gr_data.Draw("AP")
    mg.Draw("A")
    ymin = sm_xsec - 5*sm_xsec_frac_error*sm_xsec
    ymax = sm_xsec + 10*sm_xsec_frac_error*sm_xsec
    mg.SetMinimum(ymin)
    mg.SetMaximum(ymax)
    mg.GetYaxis().CenterTitle()
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetYaxis().SetTitleOffset(1.05)
    mg.GetYaxis().SetLabelSize(0.05)
    mg.GetXaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetTitleOffset(1.1)
    mg.GetXaxis().SetLabelSize(0.05)
    mg.GetXaxis().SetRangeUser(xmin+1,xmax-1)
    
    
    l = ROOT.TLegend(0.35,0.6,0.94,0.88)
    l.SetBorderSize(1)
    l.AddEntry(gr_data,"sample points","pe1")
    l.AddEntry(gr_fitted_xsec,"fitted cross section","l")
    l.AddEntry(gr_CMS,"CMS result 95% CL (2.3 fb^{-1})","fl")
    l.AddEntry(gr_CMS_300fb,"CMS prospect 95% CL (300 fb^{-1})","fl")
    
    l.Draw("same")
    
    
    downpad.cd()
    ROOT.gPad.SetMargin(0.15,0.05,0.25,0.01)
    y_1 = [(i-sm_xsec)**2/((sm_xsec_frac_error*sm_xsec)**2) for i in y]
    ymin = -2
    ymax = 15
    gr_chi2 = ROOT.TGraph(len(x),array('d',x),array('d',y_1))
    gr_chi2.SetLineWidth(2)
    gr_chi2.SetLineColor(1)
    gr_chi2.SetLineStyle(1)
    # prospects @ 300 fb
    y_2 = [(i-min(y))**2/((sm_xsec_frac_error_300fb*min(y))**2) for i in y]
    gr_chi2_300fb = ROOT.TGraph(len(x),array('d',x),array('d',y_2))
    gr_chi2_300fb.SetLineWidth(2)
    gr_chi2_300fb.SetLineColor(1)
    gr_chi2_300fb.SetLineStyle(2)
    
    
    # get roots
    p_1 = P.fit(x, y_1, 4)
    roots_1 = sorted([i.real for i in (p_1 - 3.84).roots() if i.imag == 0])
    print roots_1
    while len(roots_1)>2: roots_1 = roots_1[1:-1]
    x_limit_1 = np.arange(roots_1[0], roots_1[1]+abs(roots_1[1]-roots_1[0])/200., abs(roots_1[1]-roots_1[0])/200.)
    gr_CMS_limit = ROOT.TGraphErrors(len(x_limit_1),array('d',x_limit_1),array('d',[5]*len(x_limit_1)),array('d',np.zeros(len(x_limit_1))),array('d',[5000]*len(x_limit_1)))
    gr_CMS_limit.SetFillColorAlpha(3,0.5)
    gr_CMS_limit.SetLineWidth(2)
    
    p_2 = P.fit(x, y_2, 4)
    roots_2 = sorted([i.real for i in (p_2 - 3.84).roots() if i.imag == 0])
    print roots_2
    while len(roots_2)>2: roots_2 = roots_2[1:-1]
    x_limit_2 = np.arange(roots_2[0], roots_2[1]+abs(roots_2[1]-roots_2[0])/200., abs(roots_2[1]-roots_2[0])/200.)
    gr_CMS_limit_300fb = ROOT.TGraphErrors(len(x_limit_2),array('d',x_limit_2),array('d',[5]*len(x_limit_2)),array('d',np.zeros(len(x_limit_2))),array('d',[5000]*len(x_limit_2)))
    gr_CMS_limit_300fb.SetFillColorAlpha(5,0.7)
    gr_CMS_limit_300fb.SetLineWidth(2)
    
    
    mg2 = ROOT.TMultiGraph("mg2",";%s [TeV^{-2}];#chi^{2}"%convertToLatex(coupling))
    mg2.Add(gr_CMS_limit,"3")
    mg2.Add(gr_CMS_limit_300fb,"3")
    mg2.Add(gr_chi2,"l")
    mg2.Add(gr_chi2_300fb,"l")
    
    mg2.Draw("A")
    mg2.SetMinimum(ymin)
    mg2.SetMaximum(ymax)
    mg2.GetYaxis().CenterTitle()
    mg2.GetYaxis().SetTitleSize(0.085)
    mg2.GetYaxis().SetTitleOffset(0.75)
    mg2.GetYaxis().SetLabelSize(0.075)
    mg2.GetXaxis().CenterTitle()
    mg2.GetXaxis().SetTitleSize(0.085)
    mg2.GetXaxis().SetTitleOffset(1.3)
    mg2.GetXaxis().SetLabelSize(0.085)
    mg2.GetXaxis().SetRangeUser(xmin+1,xmax-1)
    
    # line
    line = ROOT.TLine(xmin+1,3.84,xmax-1,3.84)
    line.SetLineColor(13)
    line.SetLineWidth(1)
    line.SetLineStyle(2)
    line.Draw("same")
    
    
    limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots_1[0],roots_1[1]))
    limits_file_300fb.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    
    c.SaveAs('%s/fitted_xsec/fit_xsec_%s.png'%(args.InputDir,coupling))
    c.SaveAs('%s/fitted_xsec/fit_xsec_%s.pdf'%(args.InputDir,coupling))
    
    print "done, saved in %s/fitted_xsec/fit_xsec_%s.png"%(args.InputDir,coupling)
    
limits_file.close()
limits_file_300fb.close()    

    
    # f_, axarr_ = plt.subplots(2, sharex=True, figsize=(10,8))
#     
#     #plot fitted function
#     xmin = -40
#     xmax = +40
#     nsteps = 50
#     x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
#     #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
#     y = func(x,*(coeff))
#     y_up = func(x,*(coeff + [i for i in errors]))
#     y_down = func(x,*(coeff - [i for i in errors]))
#     axarr_[0].plot(x,y,c='r',label = 'fit')
#     axarr_[0].errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b",label='sample points')
#     #axarr_[0].fill_between(x, y_down,y_up, facecolor='green', alpha=0.3, label=r'fit $\pm 1\sigma$')
#     axarr_[0].fill_between(x, sm_xsec - sm_xsec_frac_error*sm_xsec, sm_xsec + sm_xsec_frac_error*sm_xsec, facecolor='green', alpha=0.3, label=r'SM 68\% CL')
#     axarr_[0].fill_between(x, sm_xsec - 2*sm_xsec_frac_error*sm_xsec, sm_xsec + 2*sm_xsec_frac_error*sm_xsec, facecolor='yellow', alpha=0.3, label=r'SM 95\% CL')
#     ymin = min(xsec)-(max(xsec)- min(xsec))/5.
#     ymax = sm_xsec + 3*sm_xsec_frac_error*sm_xsec
#     axarr_[0].axis([xmin, xmax, ymin, ymax])
#     axarr_[0].grid(True)
#     axarr_[0].set_xlabel(coupling, fontsize = 15)
#     axarr_[0].set_ylabel(r'cross section [pb]', fontsize = 15)
#     axarr_[0].legend(loc="upper right")
#     axarr_[0].text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=20, color="r")
#     
#     y = [(i-sm_xsec)**2/((sm_xsec_frac_error*sm_xsec)**2) for i in y]
#     ymin = -2
#     ymax = 10
#     axarr_[1].plot(x,y,c="r",label=r'fitted $\chi^{2}$')
#     axarr_[1].axhline(y=3.84, linestyle="dashed", linewidth=2, color="navy", label=r'95\% CL')
#     p = P.fit(x, y, 4)
#     roots = sorted([i.real for i in (p - 3.84).roots() if i.imag == 0])
#     print roots
#     while len(roots)>2: roots = roots[1:-1]
#     axarr_[1].axvline(x=roots[0], linestyle="dashed", linewidth=1, color="navy")
#     axarr_[1].axvline(x=roots[1], linestyle="dashed", linewidth=1, color="navy")
#     axarr_[1].text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
#     axarr_[1].text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
#     axarr_[1].axis([xmin, xmax, ymin, ymax])
#     axarr_[1].grid(True)
#     axarr_[1].set_xlabel(r'$%s$'%convertToLatex(coupling), fontsize = 15)
#     axarr_[1].set_ylabel(r'$\chi^{2}$', fontsize = 15)
#     axarr_[1].legend(loc="upper right")
    
    # plot points with error bars
    #plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
    # plt.axis([xmin, xmax, ymin, ymax])
#     plt.grid(True)
#     plt.xlabel(coupling, fontsize = 15)
#     plt.ylabel('cross section [pb]', fontsize = 15)
    
    # draw some things on the canvas
    #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
    #axarr_[0].text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=20, color="r")
    #plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(nevents[0]), fontsize=17, color="b")
    
    #print '%s: sigma=%.4f ( 1 + %.4fC + %.4f C^2 )$'%(coupling,coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0])
    # print '%s: sigma=%.4f ( 1 + %.4fC + %.4f C^2 )$'%(coupling,coeff[0],coeff[1],coeff[2])
#     
#     plt.subplots_adjust(wspace=0.02, hspace=0)
#     
#     f_.savefig('%s/fitted_xsec/fit_xsec_%s.png'%(args.InputDir,coupling))
#     f_.savefig('%s/fitted_xsec/fit_xsec_%s.pdf'%(args.InputDir,coupling))
#     
#     f_.clf()


   # print "done, saved in %s/fitted_xsec/fit_xsec_%s.png"%(args.InputDir,coupling)

# fsm_ = open(args.SMxsecDir+ "/cross_section_SM.txt", 'r')
# lines = fsm_.readlines()
# sm_xsec =  float(lines[-1].split(" ")[2])
# sm_xsec_error = float(lines[-1].split(" ")[3])
# fsm_.close()

# f, axarr = plt.subplots(3, sharex=True, figsize=(12,8))
# x_ = range(len(coupling_names))
# latex_coupling_names = [r'$%s$'%convertToLatex(c) for c in coupling_names]
# plt.xticks(x_, latex_coupling_names)
# axarr[0].errorbar(x_, p0, yerr = p0e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[0].set_ylabel("p0 [pb]", fontsize = 15)
# axarr[0].grid(True)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].set_ylim(sm_xsec - 3*sm_xsec_error, sm_xsec + 3*sm_xsec_error)
# textymin, textymax = axarr[0].get_ylim()
# axarr[0].text(x_[0]-1,textymax + (textymax - textymin)/float(10), r'$\sigma=p0 ( 1 + p1 \ C + p2 \ C^2 )$', fontsize=20)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].fill_between([-1] + x_ + [len(x_)+1], sm_xsec - (sm_xsec_frac_error*sm_xsec), sm_xsec + (sm_xsec_frac_error*sm_xsec), facecolor='red', alpha=0.3, label='SM cross section')
# axarr[0].legend(loc='lower right')
# axarr[1].errorbar(x_, p1, yerr = p1e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[1].set_ylabel("p1", fontsize = 15)
# axarr[1].grid(True)
# axarr[2].errorbar(x_, p2, yerr = p2e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[2].set_ylabel("p2", fontsize = 15)
# axarr[2].set_xlabel("non-zero EFT coupling", fontsize = 15)
# axarr[2].grid(True)
# f.savefig('%s/fitted_xsec/fit_xsec_coeff_results.png'%(args.InputDir))
# f.savefig('%s/fitted_xsec/fit_xsec_coeff_results.pdf'%(args.InputDir))