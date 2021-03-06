import ROOT
import os
from argparse import ArgumentParser
import sys
from math import sqrt
from numpy import polyfit, diag, sqrt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial as P
from numpy.random import normal, poisson


def extract_xsec(dir,c,v):
    if os.path.isfile(dir + "/" + c + "/" + c + v + "/cross_section.txt"):
        f_ = open(dir + "/" + c + "/" + c + v + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec = float(lines[-1].split(" ")[2])
        return xsec
    else:
        print "ERROR: could not extract cross section from "+args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt"
        sys.exit(1)

def get_Original_nevents(dir,c,v,nevents_dict):
    result=0
    files = [f for f in os.listdir(dir) if c+v in f]
    if len(files)==0:
        print "ERROR: couldn't retrieve original number of events for coupling %s, %s"%(c,v)
        sys.exit(1)
    for run,n in nevents_dict.iteritems():
        for f in files:
            if run in f: result += n
    
    return result


def convert_coupling_toNumerical(v):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in v: return -float(v.replace("minus","").replace("p","."))
    else: return float(v.replace("p","."))

        
def func(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)        


def fitTemplate(data_hist,EFT_templ,SM_templ,coupl_name, coupl_value, savedir, xaxis_title = "NN output",SM_unc=0.2):

    #########################
    #
    # Add uncertainty to the data
    #
    #########################
    frac_error = SM_unc 
    for bin in range(data_hist.GetNbinsX()):
        content = data_hist.GetBinContent(bin+1)
        if content != 0:
            random_var = normal(0,frac_error*content)
            data_hist.SetBinContent(bin+1,content+random_var)
            data_hist.SetBinError(bin+1,data_hist.GetBinError(bin+1)+data_hist.GetBinContent(bin+1)*frac_error)
        
    #########################
    #
    # Plot the prefitted histos/teplates
    #
    #########################
    # c_prefit = ROOT.TCanvas("c_prefit","c_prefit",600,500)
#     c_prefit.cd()
#     ROOT.gPad.SetMargin(0.15,0.1,0.15,0.1)
#     EFT_templ.SetFillColor(4)
#     EFT_templ.SetLineWidth(0)
#     #EFT_templ.Draw("hist")
#     SM_templ.SetFillColor(2)
#     SM_templ.SetLineWidth(0)
#     #SM_templ.Draw("same hist")
#     stack_prefit = ROOT.THStack("stack_prefit","")
#     stack_prefit.Add(EFT_templ)
#     stack_prefit.Add(SM_templ)
#     stack_prefit.SetMinimum(0)
#     stack_prefit.SetMaximum(2*stack_prefit.GetMaximum())
#     stack_prefit.Draw("hist")
#     stack_prefit.GetXaxis().SetTitle(xaxis_title)
#     stack_prefit.GetXaxis().SetTitleOffset(1.2)
#     stack_prefit.GetXaxis().SetTitleSize(0.055)
#     stack_prefit.GetYaxis().SetTitle("Entries")
#     stack_prefit.GetYaxis().SetTitleOffset(1.2)
#     stack_prefit.GetYaxis().SetTitleSize(0.055)
#     data_hist.SetMarkerStyle(20)
#     data_hist.SetMarkerSize(1)
#     data_hist.SetLineWidth(2)
#     data_hist.Draw("same pX1E1")
#     l_prefit = ROOT.TLegend(0.5,0.65,0.9,0.9)
#     l_prefit.SetHeader("prefit")
#     l_prefit.AddEntry(EFT_templ,"EFT template","f")
#     l_prefit.AddEntry(SM_templ,"SM template","f")
#     l_prefit.AddEntry(data_hist,"data","pe")
#     l_prefit.Draw("same")
#     c_prefit.SaveAs("prefit.pdf")
    
    
    
    #print SM_templ.GetXaxis().GetBinLowEdge(1), SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX())
    var = ROOT.RooRealVar("var","var",SM_templ.GetXaxis().GetBinLowEdge(1),SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX()))
    SM_Hist = ROOT.RooDataHist("SM_Hist","SM_Hist",ROOT.RooArgList(var),SM_templ)
    SM_PDF = ROOT.RooHistPdf("SM_PDF","SM_PDF",ROOT.RooArgSet(var),SM_Hist)
    EFT_Hist = ROOT.RooDataHist("EFT_Hist","EFT_Hist",ROOT.RooArgList(var),EFT_templ)
    EFT_PDF = ROOT.RooHistPdf("EFT_PDF","EFT_PDF",ROOT.RooArgSet(var),EFT_Hist)
    data_Hist = ROOT.RooDataHist("data_Hist","data_Hist",ROOT.RooArgList(var),data_hist)
    data_PDF = ROOT.RooHistPdf("data_PDF","data_PDF",ROOT.RooArgSet(var),data_Hist)
    
    
    
    N_SM = ROOT.RooRealVar("N_SM","N_SM",data_hist.Integral(),0,5*(data_hist.Integral())) 
    N_EFT = ROOT.RooRealVar("N_EFT","N_EFT",data_hist.Integral(),0,5*(data_hist.Integral()))
    mod = ROOT.RooAddPdf("mod","mod",ROOT.RooArgList(SM_PDF,EFT_PDF),ROOT.RooArgList(N_SM,N_EFT)); 
    
    fitRes = mod.fitTo(data_Hist,ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-3), ROOT.RooFit.Verbose(False), ROOT.RooFit.Extended(True), ROOT.RooFit.Save())

    fitted_EFT = N_EFT.getVal()
    sigma_fitted_EFT = N_EFT.getError()
    fitted_SM = N_SM.getVal()
    sigma_fitted_SM = N_SM.getError()

    print fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM
    
    #########################
    #
    # Plot the postfitted histos/teplates
    #
    #########################
    EFT_templ_copy = EFT_templ.Clone()
    SM_templ_copy = SM_templ.Clone()
    c_postfit = ROOT.TCanvas("c_postfit","c_postfit",600,500)
    c_postfit.cd()
    ROOT.gPad.SetMargin(0.15,0.1,0.15,0.1)
    EFT_templ_copy.SetFillColor(4)
    EFT_templ_copy.SetLineWidth(0)
    SM_templ_copy.SetFillColor(2)
    SM_templ_copy.SetLineWidth(0)
    EFT_templ_copy.Scale(fitted_EFT/EFT_templ_copy.Integral())
    SM_templ_copy.Scale(fitted_SM/SM_templ_copy.Integral())
    stack_postfit = ROOT.THStack("stack_postfit","")
    stack_postfit.Add(EFT_templ_copy)
    stack_postfit.Add(SM_templ_copy)
    stack_postfit.SetMinimum(0)
    stack_postfit.SetMaximum(2.5*stack_postfit.GetMaximum())
    stack_postfit.Draw("hist")
    stack_postfit.GetXaxis().SetTitle(xaxis_title)
    stack_postfit.GetXaxis().SetTitleOffset(1.2)
    stack_postfit.GetXaxis().SetTitleSize(0.055)
    stack_postfit.GetYaxis().SetTitle("Entries")
    stack_postfit.GetYaxis().SetTitleOffset(1.2)
    stack_postfit.GetYaxis().SetTitleSize(0.055)
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(1)
    data_hist.SetLineWidth(2)
    data_hist.Draw("same pX1E1")
    l_postfit = ROOT.TLegend(0.5,0.65,0.9,0.9)
    l_postfit.SetHeader("postfit")
    l_postfit.AddEntry(EFT_templ_copy,"EFT template","f")
    l_postfit.AddEntry(SM_templ_copy,"SM template","f")
    l_postfit.AddEntry(data_hist,"data","pe")
    l_postfit.Draw("same")
    t = ROOT.TPaveText(0.16,0.7,0.49,0.89,"NDC")
    t.SetTextSize(0.04)
    t.SetTextFont(42)
    t.AddText("N_{EFT} = %i #pm %i"%(int(fitted_EFT),int(sigma_fitted_EFT)))
    t.AddText("N_{SM} = %i #pm %i"%(int(fitted_SM),int(sigma_fitted_SM)))
    t.Draw("same")
    t2 = ROOT.TPaveText(0.14,0.91,0.6,0.97,"NDC")
    t2.SetTextSize(0.04)
    t2.SetTextFont(42)
    t2.AddText("coupling: %s, value: %s"%(coupl_name, coupl_value))
    t2.Draw("same")
    c_postfit.SaveAs("%s/postfit_%s_%s.pdf"%(savedir,coupl_name, coupl_value))
    #print "Saved fitted output histograms in %s/postfit_%s_%s.pdf"%(savedir,coupl_name, coupl_value)
    
    return fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM
    
    

def main():
    ROOT.gROOT.SetBatch(True)
    
    parser = ArgumentParser()
    parser.add_argument('--TemplateFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/templates/templates.root", help='root file with templates')
    parser.add_argument('--ValidationFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/validation_data.root", help='root file with validation data')
    parser.add_argument('--OutputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/control_plots/", help='path to the converted delphes training files')
    parser.add_argument('--XsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttcc_DiLepton_ValidationPerCoupling_WithInterference/", help='path to MODELSCAN directory that holds the cross sections')
    parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/", help='path to MODELSCAN directory that holds the cross sections')
    args = parser.parse_args()
    
    if not os.path.isdir(args.OutputDir): os.mkdir(args.OutputDir)
    
    f_templ = ROOT.TFile(args.TemplateFile)
    f_data = ROOT.TFile(args.ValidationFile)
    
    fitting_dict = { #name:class_number
            #'SM': {},
            #LLLL
            'C11Qq': {},
            'C81Qq': {},
            'C13Qq': {},
            'C83Qq': {},
            #RRRR
            'C1tu': {},
            'C8tu': {},
            'C1td': {},
            'C8td': {},
             #LLRR
            'C8tq': {},
            'C8Qu': {},
            'C8Qd': {},
            'C1tq': {},
            'C1Qu': {},
            'C1Qd': {}    
        }
    
    fracerr_measured_sm=0.2
    limits_file = open('%s/limits_templateFit.txt'%(args.XsecDir),"w")
    #limits_file = open('%s/limits_nocut.txt'%(args.InputDir),"w")
    limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))

    
    for c,result_dict in fitting_dict.iteritems():
        
        coupling_handedness=""
        if "Qq" in c: coupling_handedness="LLLL"
        elif "Qu" in c or "Qd" in c or "tq" in c: coupling_handedness="LLRR"
        else: coupling_handedness="RRRR"
        
        # loop over histograms in validation file
        histograms = []
        for key in f_data.GetListOfKeys():
            hist = key.ReadObj()
            if c in hist.GetName() and coupling_handedness in hist.GetName(): histograms.append(hist)  
        
        #extract templates
        templates = {}
        for key in f_templ.GetListOfKeys():
            hist = key.ReadObj()
            if coupling_handedness in hist.GetName() and "SM" in hist.GetName(): templates["SM"] = hist
            if c in hist.GetName() and coupling_handedness in hist.GetName(): 
                v = hist.GetName().split("_")[2]
                templates[c+v] = hist
        
        
        xsec = {} #in fb
        original_nevents = {}
        nevents={"run_01":10000,"run_02":20000}
        
        #First do the SM (and extract xsections)
        for h in histograms:
            name = h.GetName()
            coupl_value = name.split("_")[2]
            xsec[coupl_value] = extract_xsec(args.XsecDir,c,coupl_value)*1000
            original_nevents[coupl_value] = get_Original_nevents(args.ValidationDir,c,coupl_value,nevents)
            if coupl_value == "0p0":
                fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM = fitTemplate(h,templates[c+"1p0"],templates["SM"],c,coupl_value,args.OutputDir)
                fitting_dict[c][coupl_value]=[fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM,xsec[coupl_value],original_nevents[coupl_value]]
        
        
        # now scan the coupling values
        for h in histograms:
            name = h.GetName()
            coupl_value = name.split("_")[2]
            if coupl_value == "0p0":continue
            fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM = fitTemplate(h,templates[c+coupl_value],templates["SM"],c,coupl_value,args.OutputDir)
            fitting_dict[c][coupl_value]=[fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM,xsec[coupl_value],original_nevents[coupl_value]]
        
        int_lumi = 100 #fb-1
        SM_expected = fitting_dict[c]["0p0"][2]*int_lumi*fitting_dict[c]["0p0"][4]/fitting_dict[c]["0p0"][5]
        observed = {}
        observed_unc = {}
        chi2 = {}
        couplings_numerical = {}
        for v,results in fitting_dict[c].iteritems():
            #if v == "0p0": continue
            observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][0])*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5]
            #observed[v] = (fitting_dict[c][v][0])*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5]
            observed_unc[v] = sqrt((fitting_dict[c][v][1]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])**2 + (fitting_dict[c][v][3]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])**2)
            #observed_unc[v] = (fitting_dict[c][v][1]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])
            #chi2[v] = (observed[v])**2/observed_unc[v]**2
            #chi2[v] = (observed[v]-SM_expected)**2/observed_unc[v]**2
            chi2[v] = (observed[v]-SM_expected)**2/SM_expected
            couplings_numerical[v]=convert_coupling_toNumerical(v)
        
        chi2_array = []
        c_array= []
        chi2Error_array=[]
        for v,res in chi2.iteritems():
               chi2_array.append(res) 
               c_array.append(couplings_numerical[v])
               chi2Error_array.append((2*abs(observed[v]-SM_expected)/SM_expected)*observed_unc[v])
        
        print ""
        print SM_expected
        print observed
        print observed_unc
        print chi2
        print ""
        
        plt.errorbar(c_array, chi2_array, yerr=chi2Error_array, fmt='o')
        
        
        # start fitting the chi2
        coeff, covmat = curve_fit(func,c_array,chi2_array,sigma=chi2Error_array)#,sigma=[0]*len(c_array))#absolute_sigma
        errors = sqrt(diag(covmat))
        
        xmin = -12
        xmax = +12
        ymin = 0
        ymax = max(chi2_array)
        nsteps = 30
        x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
        #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
        y = coeff[0]+ coeff[1]*x + coeff[2]*x*x + coeff[3]*x*x*x + + coeff[4]*x*x*x*x
        #yup = (coeff[0]+errors[0])*(1 + (coeff[1]+errors[1])*x + (coeff[2]+errors[2])*x*x)
        #ydown = (coeff[0]-errors[0])*(1 + (coeff[1]-errors[1])*x + (coeff[2]-errors[2])*x*x)
        
        plt.plot(x,y,c="r",label='fit')
        plt.axhline(y=3.841, linestyle="dashed", linewidth=2, color="navy")
        p = P.fit(x, y, 4)
        roots = sorted([i.real for i in (p - 3.841).roots() if i.imag == 0])
        print roots
        while len(roots)>2: roots = roots[1:-1]
        #plt.axvline(x=roots[0], linestyle="dashed", linewidth=2, color="navy")
        #plt.axvline(x=roots[1], linestyle="dashed", linewidth=2, color="navy")
        plt.legend(loc="upper right")
    
    
        # plot points with error bars
        # plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
        #if roots[0] < xmin: plt.axis([roots[0]-1, roots[1]+1, ymin, ymax])
        #else: 
        plt.axis([xmin, xmax, ymin, ymax])
        plt.grid(True)
        plt.xlabel(c, fontsize = 15)
        plt.ylabel('#chi^{2} value', fontsize = 15)
    
        # draw some things on the canvas
        #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
        #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.5f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=12, color="r")
        #plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(sum(nevents.values())), fontsize=17, color="b")
        #plt.text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
        #plt.text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
        
        plt.savefig('%s/fit_Chi2_%s.png'%(args.OutputDir,c))
        plt.savefig('%s/fit_Chi2_%s.pdf'%(args.OutputDir,c))
        
        plt.cla()
        
        #limits_file.write("%s & [%.2f,%.2f] \n"%(c,roots[0],roots[1]))  
        
        
    limits_file.close()    
        
    #print fitting_dict

    
    #fitTemplate(inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root", inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order2_run_01_tag_1_delphes_events.root",inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root",template_hist_name="hist_top_eta")
    
    # f_data = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/validation_data.root")
#     data_hist = ROOT.TH1D()
#     data_hist = f_data.Get("h_C1tq_10p0_outLLRR")
#     f_EFT_templ = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/templates/templates.root")
#     EFT_templ = ROOT.TH1D()
#     EFT_templ = f_EFT_templ.Get("h_C1tq_outLLRR")
#     f_SM_templ = ROOT.TFile("/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Training_WithInterference/templates/templates.root")
#     SM_templ = ROOT.TH1D()
#     SM_templ = f_SM_templ.Get("h_SM_outLLRR")
#     #SM_templ.Scale(10.)
#     
#     outdir="/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ValidationPerCoupling_WithInterference/templates/"
#     fitTemplate(data_hist,EFT_templ,SM_templ,"C1tq","10p0",outdir)


if __name__ == "__main__":
    main()