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
from numpy.random import normal


def extract_xsec(dir,c,v):
    if os.path.isfile(dir + "/" + c + "/" + c + "_" + v + "/cross_section.txt"):
        f_ = open(dir + "/" + c + "/" + c + "_" + v + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec = float(lines[-1].split(" ")[2])
        return xsec
    else:
        print "ERROR: could not extract cross section from "+args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt"
        sys.exit(1)

def get_Original_nevents(dir,c,v,nevents_dict):
    result=0
    files = [f for f in os.listdir(dir) if c+"_"+v in f]
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

        
def func2(x, a, b, c):
        return (a + b*x + c*x*x)   

def func4(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)      


def fitTemplate(data_hist,EFT_templ,SM_templ,coupl_name, coupl_value, savedir, xaxis_title = "NN output",SM_unc=0.2,npseudo=100):
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    #ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 2001;")
    #m = ROOT.RooMinuit()
    #m.setPrintLevel(-1)
    #########################
    #
    # Add uncertainty to the data
    #
    #########################
    fitted_EFT_arr = []
    fitted_SM_arr = []
    fitted_chi2_array = []
    for i in range(npseudo):
        pseudo_data_hist = data_hist.Clone()
        frac_error = SM_unc 
        for bin in range(pseudo_data_hist.GetNbinsX()):
            content = pseudo_data_hist.GetBinContent(bin+1)
            if content != 0:
                random_var = normal(0,frac_error*content)
                pseudo_data_hist.SetBinContent(bin+1,content+random_var)
                pseudo_data_hist.SetBinError(bin+1,pseudo_data_hist.GetBinError(bin+1)+pseudo_data_hist.GetBinContent(bin+1)*frac_error)
    
    
    
        #print SM_templ.GetXaxis().GetBinLowEdge(1), SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX())
        var = ROOT.RooRealVar("var","var",SM_templ.GetXaxis().GetBinLowEdge(1),SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX()))
        SM_Hist = ROOT.RooDataHist("SM_Hist","SM_Hist",ROOT.RooArgList(var),SM_templ)
        SM_PDF = ROOT.RooHistPdf("SM_PDF","SM_PDF",ROOT.RooArgSet(var),SM_Hist)
        EFT_Hist = ROOT.RooDataHist("EFT_Hist","EFT_Hist",ROOT.RooArgList(var),EFT_templ)
        EFT_PDF = ROOT.RooHistPdf("EFT_PDF","EFT_PDF",ROOT.RooArgSet(var),EFT_Hist)
        data_Hist = ROOT.RooDataHist("pseudo_data_hist","pseudo_data_hist",ROOT.RooArgList(var),pseudo_data_hist)
        data_PDF = ROOT.RooHistPdf("data_PDF","data_PDF",ROOT.RooArgSet(var),data_Hist)
    
    
    
        N_SM = ROOT.RooRealVar("N_SM","N_SM",pseudo_data_hist.Integral(),0,5*(pseudo_data_hist.Integral())) 
        N_EFT = ROOT.RooRealVar("N_EFT","N_EFT",pseudo_data_hist.Integral(),0,5*(pseudo_data_hist.Integral()))
        mod = ROOT.RooAddPdf("mod","mod",ROOT.RooArgList(SM_PDF,EFT_PDF),ROOT.RooArgList(N_SM,N_EFT)); 
    
        fitRes = mod.fitTo(data_Hist,ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1000), ROOT.RooFit.Verbose(False), ROOT.RooFit.Extended(True), ROOT.RooFit.Save())

        fitted_EFT = N_EFT.getVal()
        sigma_fitted_EFT = N_EFT.getError()
        fitted_SM = N_SM.getVal()
        sigma_fitted_SM = N_SM.getError()
        
        fitted_EFT_arr.append(fitted_EFT)
        fitted_SM_arr.append(fitted_SM)
        
        #Chi2 fit
        fitted_chi2 = ROOT.RooChi2Var("fitted_chi2","fitted_chi2",mod,data_Hist)#,ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1000), ROOT.RooFit.Verbose(False), ROOT.RooFit.Extended(False), ROOT.RooFit.Save())
        fitted_chi2_array.append(fitted_chi2.getVal())
    
    mean_fitted_EFT = np.mean(np.asarray(fitted_EFT_arr))
    err_fitted_EFT = np.std(np.asarray(fitted_EFT_arr))
    mean_fitted_SM = np.mean(np.asarray(fitted_SM_arr))
    err_fitted_SM = np.std(np.asarray(fitted_SM_arr))
    mean_fitted_chi2 = np.mean(np.asarray(fitted_chi2_array))
    err_fitted_chi2 = np.std(np.asarray(fitted_chi2_array))
    print mean_fitted_EFT,err_fitted_EFT , mean_fitted_SM, err_fitted_SM, mean_fitted_chi2, err_fitted_chi2
    
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
    pseudo_data_hist.SetMarkerStyle(20)
    pseudo_data_hist.SetMarkerSize(1)
    pseudo_data_hist.SetLineWidth(2)
    pseudo_data_hist.Draw("same pX1E1")
    l_postfit = ROOT.TLegend(0.5,0.65,0.9,0.9)
    l_postfit.SetHeader("postfit")
    l_postfit.AddEntry(EFT_templ_copy,"EFT template","f")
    l_postfit.AddEntry(SM_templ_copy,"SM template","f")
    l_postfit.AddEntry(pseudo_data_hist,"data","pe")
    l_postfit.Draw("same")
    t = ROOT.TPaveText(0.16,0.7,0.49,0.89,"NDC")
    t.SetTextSize(0.04)
    t.SetTextFont(42)
    t.AddText("N_{EFT} = %i #pm %i"%(int(fitted_EFT),int(sigma_fitted_EFT)))
    t.AddText("N_{SM} = %i #pm %i"%(int(fitted_SM),int(sigma_fitted_SM)))
    t.AddText("#chi^{2}/ndof = %i #pm %i"%(mean_fitted_chi2,err_fitted_chi2))
    t.Draw("same")
    t2 = ROOT.TPaveText(0.14,0.91,0.6,0.97,"NDC")
    t2.SetTextSize(0.04)
    t2.SetTextFont(42)
    t2.AddText("coupling: %s, value: %s"%(coupl_name, coupl_value))
    t2.Draw("same")
    c_postfit.SaveAs("%s/postfit_%s_%s.pdf"%(savedir,coupl_name, coupl_value))
    print "Saved fitted output histograms in %s/postfit_%s_%s.pdf"%(savedir,coupl_name, coupl_value)
    
    return mean_fitted_EFT,err_fitted_EFT , mean_fitted_SM, err_fitted_SM, mean_fitted_chi2, err_fitted_chi2
    
    

def main():
    ROOT.gROOT.SetBatch(True)
    
    parser = ArgumentParser()
    parser.add_argument('--TemplateFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_TRAINING/templates/templates.root", help='root file with templates')
    parser.add_argument('--ValidationFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings/templates/validation_data.root", help='root file with validation data')
    parser.add_argument('--OutputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings/templates/control_plots/", help='path to the converted delphes training files')
    parser.add_argument('--XsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings/", help='path to MODELSCAN directory that holds the cross sections')
    parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings/", help='path to MODELSCAN directory that holds the cross sections')
    args = parser.parse_args()
    
    if not os.path.isdir(args.OutputDir): os.mkdir(args.OutputDir)
    
    f_templ = ROOT.TFile(args.TemplateFile)
    f_data = ROOT.TFile(args.ValidationFile)
    
    fitting_dict = { #name:class_number
            #'SM':  {},
            #LLLL
            #"cQQ1": {},
            "cQQ8": {},
            #LLRR
            "cQt1": {},
            "cQb1": {},
            "cQt8": {},
            "cQb8": {},
            #RRRR
            "ctb1": {},
            "ctb8": {},
        }
    
    fracerr_measured_sm=0.2
    limits_file = open('%s/limits_ttbb_templateFit_chi2.txt'%(args.XsecDir),"w")
    #limits_file = open('%s/limits_nocut.txt'%(args.InputDir),"w")
    limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))

    
    for c,result_dict in fitting_dict.iteritems():
        if c == "cQQ1": continue
        print "Processing: %s"%c
        
        coupling_handedness=""
        if "QQ" in c: coupling_handedness="LLLL"
        elif "Qt" in c or "Qb" in c: coupling_handedness="LLRR"
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
        nevents={"run_01":30000,"run_02":0}
        
        #First do the SM (and extract xsections)
        for h in histograms:
            name = h.GetName()
            coupl_value = name.split("_")[2]
            print coupl_value
            #print args.XsecDir,c,coupl_value
            xsec[coupl_value] = extract_xsec(args.XsecDir,c,coupl_value)*1000
            original_nevents[coupl_value] = get_Original_nevents(args.ValidationDir,c,coupl_value,nevents)
            if coupl_value == "0p0":
                fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM, fitted_chi2, sigma_fitted_chi2 = fitTemplate(h.Clone(),templates[c+"1p0"].Clone(),templates["SM"].Clone(),c,coupl_value,args.OutputDir)
                fitting_dict[c][coupl_value]=[fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM,xsec[coupl_value],original_nevents[coupl_value], fitted_chi2, sigma_fitted_chi2]
        
        
        # now scan the coupling values
        for h in histograms:
            name = h.GetName()
            coupl_value = name.split("_")[2]
            print "Processing %s, %s"%(name,coupl_value)
            if coupl_value == "0p0":continue
            fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM, fitted_chi2, sigma_fitted_chi2 = fitTemplate(h.Clone(),templates[c+coupl_value].Clone(),templates["SM"].Clone(),c,coupl_value,args.OutputDir)
            fitting_dict[c][coupl_value]=[fitted_EFT,sigma_fitted_EFT , fitted_SM, sigma_fitted_SM,xsec[coupl_value],original_nevents[coupl_value], fitted_chi2, sigma_fitted_chi2]
        
        int_lumi = 100 #fb-1
        SM_expected = fitting_dict[c]["0p0"][0]*int_lumi*fitting_dict[c]["0p0"][4]/fitting_dict[c]["0p0"][5]
        err_SM_expected = fitting_dict[c]["0p0"][1]*int_lumi*fitting_dict[c]["0p0"][4]/fitting_dict[c]["0p0"][5]
        print SM_expected, err_SM_expected
        
        observed = {}
        observed_unc = {}
        chi2 = {}
        err_chi2 = {}
        couplings_numerical = {}
        for v,results in fitting_dict[c].iteritems():
            if v == "0p0": continue
            #observed[v] = (fitting_dict[c][v][2]+fitting_dict[c][v][0])*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5]
            observed[v] = (fitting_dict[c][v][0])*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5]
            #observed_unc[v] = sqrt((fitting_dict[c][v][1]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])**2 + (fitting_dict[c][v][3]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])**2)
            observed_unc[v] = (fitting_dict[c][v][1]*int_lumi*fitting_dict[c][v][4]/fitting_dict[c][v][5])
            #chi2[v] = fitting_dict[c][v][6]
            err_chi2[v] = fitting_dict[c][v][7]
            chi2[v] = (observed[v]-SM_expected)**2/observed_unc[v]**2
            couplings_numerical[v]=convert_coupling_toNumerical(v)
        
        obs_low95CL = []
        obs_high95CL = []
        observed_arr = []
        unc_arr = []
        c_array= []
        chi2_arr = []
        err_chi2_arr = []
        for v,res in observed.iteritems():
            observed_arr.append(res)
            obs_low95CL.append(res-2*observed_unc[v]) 
            obs_high95CL.append(res+2*observed_unc[v])
            c_array.append(couplings_numerical[v]) 
            unc_arr.append(observed_unc[v])
            chi2_arr.append(chi2[v])
            err_chi2_arr.append(err_chi2[v])
        
        #plt.errorbar(c_array, chi2_arr, yerr=[0]*len(chi2_arr), fmt='o', label='samples data points')
        
        coeff_centr, covmat_centr = curve_fit(func2,c_array,observed_arr,sigma=unc_arr)#,sigma=[0]*len(c_array))#absolute_sigma
        errors_centr = sqrt(diag(covmat_centr))
        #coeff_centr[0] = coeff_centr[0]-SM_expected
        
        # coeff_up, covmat_up = curve_fit(func4,c_array,obs_high95CL)#,sigma=[0]*len(c_array))#absolute_sigma
#         errors_up = sqrt(diag(covmat_up))
#         coeff_down, covmat_down = curve_fit(func4,c_array,obs_low95CL)#,sigma=[0]*len(c_array))#absolute_sigma
#         errors_down = sqrt(diag(covmat_down))
        
        xmin = -10
        xmax = +10
       # ymax = 120)
        #ymin = 0#-ymax/4.
        nsteps = 50
        x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
        #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
        #y_centr = coeff_centr[0]+ coeff_centr[1]*x + coeff_centr[2]*x*x
        y_centr = func2(x,*(coeff_centr))
        y_centr = [(i-min(y_centr))**2/(err_SM_expected**2) for i in y_centr]
        ymax = 10#max(y_centr)
        ymin = 0
        #y_up = func4(x,*(coeff_up))
        #y_down = func4(x,*(coeff_down))
        #y_up = (coeff_centr[0]+errors_centr[0])*(1 + (coeff_centr[1]+errors_centr[1])*x + (coeff_centr[2]+errors_centr[2])*x*x)
        #y_down = (coeff_centr[0]-errors_centr[0])*(1 + (coeff_centr[1]-errors_centr[1])*x + (coeff_centr[2]-errors_centr[2])*x*x)
        #y_up = func4(x,*(coeff_centr + [i for i in errors_centr]))
        #y_down = func4(x,*(coeff_centr - [i for i in errors_centr]))
        
        plt.plot(x,y_centr,c="r",label='fit',linestyle="dotted")
        #plt.show()
        #plt.fill_between(x, y_down,y_up, facecolor='green', alpha=0.3, label=r'fitted 95% CL')
        plt.axhline(y=3.84, linestyle="dashed", linewidth=2, color="navy")
        p = P.fit(x, y_centr, 4)
        roots = sorted([i.real for i in (p - 3.84).roots() if i.imag == 0])
        print roots
        while len(roots)>2: roots = roots[1:-1]
        plt.axvline(x=roots[0], linestyle="dashed", linewidth=2, color="navy")
        plt.axvline(x=roots[1], linestyle="dashed", linewidth=2, color="navy")
        plt.legend(loc="upper right")
    
    
        # plot points with error bars
        #plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
        #if roots[0] < xmin: plt.axis([roots[0]-1, roots[1]+1, ymin, ymax])
        #else: 
        plt.axis([xmin, xmax, ymin, ymax])
        plt.grid(True)
        plt.xlabel(c, fontsize = 15)
        plt.ylabel('#chi^{2}/ndof', fontsize = 15)
    
        # draw some things on the canvas
        #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
        #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.5f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=12, color="r")
        #plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(sum(nevents.values())), fontsize=17, color="b")
        plt.text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
        plt.text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
        
        plt.savefig('%s/fit_ConfidenceLevelInterval_%s.png'%(args.OutputDir,c))
        plt.savefig('%s/fit_ConfidenceLevelInterval_%s.pdf'%(args.OutputDir,c))
        
        plt.cla()
        
        limits_file.write("%s & [%.2f,%.2f] \n"%(c,roots[0],roots[1]))  
       #  print ""
#         print SM_expected
#         print observed
#         print observed_unc
#         print chi2
#         print ""
        
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