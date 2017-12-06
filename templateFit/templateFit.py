import ROOT
import os

def fitTemplate(data_file,signal_templ_file,background_templ_file,template_hist_name="NN_output_singlet",xaxis_title = "NN output",SM_unc=0.2):
    #########################
    #
    # Retrieve the templates and the data
    #
    #########################
    f_data = ROOT.TFile(data_file)
    data_hist = ROOT.TH1D()
    data_hist = f_data.Get(template_hist_name)
    f_EFT_templ = ROOT.TFile(signal_templ_file)
    EFT_templ = ROOT.TH1D()
    EFT_templ = f_EFT_templ.Get(template_hist_name)
    f_SM_templ = ROOT.TFile(background_templ_file)
    SM_templ = ROOT.TH1D()
    SM_templ = f_SM_templ.Get(template_hist_name)
    SM_templ.Scale(10.)
    
    pseudo_data_hist = data_hist.Clone()
    pseudo_data_hist.Reset()
    for i in range(int(data_hist.GetEntries())):
        pseudo_data_hist.Fill(data_hist.GetRandom())
    #########################
    #
    # Add uncertainty to the data
    #
    #########################
    frac_error = SM_unc 
    for bin in range(data_hist.GetNbinsX()):
        pseudo_data_hist.SetBinError(bin,pseudo_data_hist.GetBinError(bin)+pseudo_data_hist.GetBinContent(bin)*frac_error)
        
    #########################
    #
    # Plot the prefitted histos/teplates
    #
    #########################
    c_prefit = ROOT.TCanvas("c_prefit","c_prefit",600,500)
    c_prefit.cd()
    ROOT.gPad.SetMargin(0.15,0.1,0.15,0.1)
    EFT_templ.SetFillColor(4)
    EFT_templ.SetLineWidth(0)
    #EFT_templ.Draw("hist")
    SM_templ.SetFillColor(2)
    SM_templ.SetLineWidth(0)
    #SM_templ.Draw("same hist")
    stack_prefit = ROOT.THStack("stack_prefit","")
    stack_prefit.Add(EFT_templ)
    stack_prefit.Add(SM_templ)
    stack_prefit.SetMinimum(0)
    stack_prefit.SetMaximum(2*stack_prefit.GetMaximum())
    stack_prefit.Draw("hist")
    stack_prefit.GetXaxis().SetTitle(xaxis_title)
    stack_prefit.GetXaxis().SetTitleOffset(1.2)
    stack_prefit.GetXaxis().SetTitleSize(0.055)
    stack_prefit.GetYaxis().SetTitle("Entries")
    stack_prefit.GetYaxis().SetTitleOffset(1.2)
    stack_prefit.GetYaxis().SetTitleSize(0.055)
    pseudo_data_hist.SetMarkerStyle(20)
    pseudo_data_hist.SetMarkerSize(1)
    pseudo_data_hist.SetLineWidth(2)
    pseudo_data_hist.Draw("same pX1E1")
    l_prefit = ROOT.TLegend(0.5,0.65,0.9,0.9)
    l_prefit.SetHeader("prefit")
    l_prefit.AddEntry(EFT_templ,"EFT template","f")
    l_prefit.AddEntry(SM_templ,"SM template","f")
    l_prefit.AddEntry(pseudo_data_hist,"data","pe")
    l_prefit.Draw("same")
    c_prefit.SaveAs("prefit.pdf")
    
    
    
    #print SM_templ.GetXaxis().GetBinLowEdge(1), SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX())
    var = ROOT.RooRealVar("var",template_hist_name,SM_templ.GetXaxis().GetBinLowEdge(1),SM_templ.GetXaxis().GetBinUpEdge(SM_templ.GetNbinsX()))
    SM_Hist = ROOT.RooDataHist("SM_Hist","SM_Hist",ROOT.RooArgList(var),SM_templ)
    SM_PDF = ROOT.RooHistPdf("SM_PDF","SM_PDF",ROOT.RooArgSet(var),SM_Hist)
    EFT_Hist = ROOT.RooDataHist("EFT_Hist","EFT_Hist",ROOT.RooArgList(var),EFT_templ)
    EFT_PDF = ROOT.RooHistPdf("EFT_PDF","EFT_PDF",ROOT.RooArgSet(var),EFT_Hist)
    data_Hist = ROOT.RooDataHist("data_Hist","data_Hist",ROOT.RooArgList(var),pseudo_data_hist)
    data_PDF = ROOT.RooHistPdf("data_PDF","data_PDF",ROOT.RooArgSet(var),data_Hist)
    
    
    
    N_SM = ROOT.RooRealVar("N_SM","N_SM",SM_templ.Integral(),0,5*(EFT_templ.Integral())+(SM_templ.Integral())) 
    N_EFT = ROOT.RooRealVar("N_EFT","N_EFT",EFT_templ.Integral(),0,5*(EFT_templ.Integral())+(SM_templ.Integral()))
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
    c_postfit = ROOT.TCanvas("c_postfit","c_postfit",600,500)
    c_postfit.cd()
    ROOT.gPad.SetMargin(0.15,0.1,0.15,0.1)
    EFT_templ.Scale(fitted_EFT/EFT_templ.Integral())
    SM_templ.Scale(fitted_SM/SM_templ.Integral())
    stack_postfit = ROOT.THStack("stack_postfit","")
    stack_postfit.Add(EFT_templ)
    stack_postfit.Add(SM_templ)
    stack_postfit.SetMinimum(0)
    stack_postfit.SetMaximum(2*stack_postfit.GetMaximum())
    stack_postfit.Draw("hist")
    stack_postfit.GetXaxis().SetTitle(xaxis_title)
    stack_postfit.GetXaxis().SetTitleOffset(1.2)
    stack_postfit.GetXaxis().SetTitleSize(0.055)
    stack_postfit.GetYaxis().SetTitle("Entries")
    stack_postfit.GetYaxis().SetTitleOffset(1.2)
    stack_postfit.GetYaxis().SetTitleSize(0.055)
    pseudo_data_hist.Draw("same pX1E1")
    l_postfit = ROOT.TLegend(0.5,0.65,0.9,0.9)
    l_postfit.SetHeader("postfit")
    l_postfit.AddEntry(EFT_templ,"EFT template","f")
    l_postfit.AddEntry(SM_templ,"SM template","f")
    l_postfit.AddEntry(pseudo_data_hist,"data","pe")
    l_postfit.Draw("same")
    c_postfit.SaveAs("postfit.pdf")
    
    

def main():
    ROOT.gROOT.SetBatch(True)
    inputdir = os.getcwd()
    #fitTemplate(inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root", inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order2_run_01_tag_1_delphes_events.root",inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root",template_hist_name="hist_top_eta")
    fitTemplate(inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/mixed.root", inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order2_run_01_tag_1_delphes_events.root",inputdir+"/CONVERTED_DELPHES_TEST_TOP_ETA/C1Qu1p0_Order0_run_01_tag_1_delphes_events.root",template_hist_name="hist_bl_DeltaR", xaxis_title = "#DeltaR(l,b)")


if __name__ == "__main__":
    main()