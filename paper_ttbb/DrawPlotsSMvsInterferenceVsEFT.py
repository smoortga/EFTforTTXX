from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TLegend, gROOT, gStyle, gPad
import sys
import ROOT
import os
import time
from argparse import ArgumentParser
from array import array
from math import *
import numpy as np

gROOT.SetBatch(1)
gStyle.SetOptStat(0)




axis_dict = {
    "deltaR_c1c2":"#Delta R(1^{st} add. b-jet,2^{nd} add. b-jet)",
    "deltaR_b1l1":"#Delta R(1^{st} lepton,1^{st} b-jet)",
    "deltaR_b1b2":"#Delta R(1^{st} b-jet,2^{nd} b-jet)",
    "deltaR_b2l2":"#Delta R(2^{nd} lepton,2^{nd} b-jet)",
    "deltaR_l1l2":"#Delta R(1^{st} lepton,2^{nd} lepton)",
    "m_l1l2":"m_{inv}(1^{st} lepton,2^{nd} lepton) [GeV]",
    "m_b2l2":"m_{inv}(2^{nd} lepton,2^{nd} b-jet) [GeV]",
    "m_c1c2":"m_{inv}(1^{st} add. b-jet,2^{nd} add. b-jet) [GeV]",
    "m_b1l1":"m_{inv}(1^{st} lepton,1^{st} b-jet) [GeV]",
    "m_b1b2":"m_{inv}(1^{st} b-jet,2^{nd} b-jet) [GeV]",
    "m_c1c2b1b2":"m_{inv}(jets) [GeV]",
    "m_c1c2b1b2l1l2":"m_{inv}(jets and leptons) [GeV]",
    "pT_l2":"p_{T}(2^{nd} lepton) [GeV]",
    "pT_b1":"p_{T}(1^{st} b-jet) [GeV]",
    "pT_c1":"p_{T}(1^{st} add. b-jet) [GeV]",
    "pT_c2":"p_{T}(2^{nd} add. b-jet) [GeV]",
    "pT_b2":"p_{T}(2^{nd} b-jet) [GeV]",
    "pT_l1":"p_{T}(1^{st} lepton) [GeV]"
}

parser = ArgumentParser()
#parser.add_argument('--ncpu', type=int, default=-1,help='number of CPU to use in parallel')
parser.add_argument('--xsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the converted delphes files')
parser.add_argument('--InputDir', default = "CONVERTED_DELPHES_ttbb_ForValidationPerCouplings_fixed", help='path to the converted delphes files')
parser.add_argument('--coupling', default="cQb8",help='name of coupling')
args = parser.parse_args()

def getxsec(filename):
    coupling = filename.split("_")[0]
    coupl_dir = filename.split("_")[0] + "_" +  filename.split("_")[1]
    if os.path.isfile(args.xsecDir + "/" + coupling + "/" + coupl_dir + "/cross_section.txt"):
        f_ = open(args.xsecDir + "/" + coupling + "/" + coupl_dir + "/cross_section.txt", 'r')
        lines = f_.readlines()
        xsec = float(lines[-1].split(" ")[2])
        return xsec


files = [i for i in os.listdir(args.InputDir) if ".root" in i and "8_" in i]


tree_0p0 = TChain("tree")
tree_minus1p0 = TChain("tree")
tree_1p0 = TChain("tree")
xsec_0p0 = 0
xsec_minus1p0 = 0
xsec_1p0 = 0
for f in files:
    #if not args.coupling in f: continue
    if "_0p0" in f: 
        tree_0p0.Add(os.getcwd() + "/" + args.InputDir  +"/" + f)
        xsec_0p0 = getxsec(f)*100
        print xsec_0p0
    elif "_minus10p0" in f : 
        tree_minus1p0.Add(os.getcwd() + "/" + args.InputDir  +"/" + f)
        xsec_minus1p0= getxsec(f)*100
        print xsec_minus1p0
    elif "_10p0" in f : 
        tree_1p0.Add(os.getcwd() + "/" + args.InputDir  +"/" + f)
        xsec_1p0= getxsec(f)*100
        print xsec_1p0

# file_0p0 = TFile(os.getcwd() + "/" + args.InputDir  +"/" + args.coupling + "_0p0_run_01_tag_1_delphes_events.root")
# file_minus1p0 = TFile(os.getcwd() + "/" + args.InputDir  +"/" + args.coupling + "_minus1p0_run_01_tag_1_delphes_events.root")
# file_1p0 = TFile(os.getcwd() + "/" + args.InputDir  +"/" + args.coupling + "_1p0_run_01_tag_1_delphes_events.root")
# 
# tree_0p0 = file_0p0.Get("tree")
# tree_minus1p0 = file_minus1p0.Get("tree")
# tree_1p0 =file_1p0.Get("tree")

branchnames = [i.GetName() for i in tree_0p0.GetListOfBranches()]

for b in branchnames:
        xmin = min([i.GetMinimum(b) for i in [tree_0p0,tree_minus1p0,tree_1p0]])
        xmax = max([i.GetMaximum(b) for i in [tree_0p0,tree_minus1p0,tree_1p0]])
        nbins = 20
        
        # SM
        sm_hist = TH1D("sm_hist",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_0p0.Draw("%s >> sm_hist"%b)
        
        #interference 
        int_hist = TH1D("int_hist",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        h_minus1p0 = TH1D("h_minus1p0",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_minus1p0.Draw("%s >> h_minus1p0"%b)
        h_1p0 = TH1D("h_1p0",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_1p0.Draw("%s >> h_1p0"%b)
        h_1p0.Scale(xsec_1p0)
        h_1p0.Add(h_minus1p0,-1*xsec_minus1p0)
        int_hist = TH1D(h_1p0)
        
        # pure EFT 
        eft_hist = TH1D("eft_hist",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        h_minus1p0_2 = TH1D("h_minus1p0_2",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_minus1p0.Draw("%s >> h_minus1p0_2"%b)
        h_1p0_2 = TH1D("h_1p0_2",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_1p0.Draw("%s >> h_1p0_2"%b)
        h_0p0_2 = TH1D("h_0p0_2",";%s;a.u."%axis_dict[b],nbins,xmin,xmax)
        tree_0p0.Draw("%s >> h_0p0_2"%b)
        h_1p0_2.Scale(xsec_1p0)
        h_1p0_2.Add(h_minus1p0_2,xsec_minus1p0)
        h_1p0_2.Add(h_0p0_2,-2*xsec_0p0)
        eft_hist = TH1D(h_1p0_2)
        
        sm_hist.SetLineColor(2)
        sm_hist.SetLineWidth(2)
        sm_hist.Scale(1./sm_hist.Integral())
        eft_hist.SetLineColor(4)
        eft_hist.SetLineWidth(2)
        eft_hist.Scale(1./eft_hist.Integral())
        int_hist.SetLineColor(8)
        int_hist.SetLineWidth(2)
        int_hist.Scale(1./int_hist.Integral())
     
        c = TCanvas("c","c",800,500)
        c.SetMargin(0.15,0.1,0.15,0.1)
        gPad.SetLogy(1)
        eft_hist.Draw("hist")
        sm_hist.Draw("hist same")
        int_hist.Draw("hist same")
        
        #style
        eft_hist.GetYaxis().SetRangeUser(0.0001,1)
        eft_hist.GetYaxis().SetTickSize(0)
        eft_hist.GetYaxis().SetLabelSize(0)
        eft_hist.GetYaxis().SetTitleSize(0.06)
        eft_hist.GetYaxis().SetTitleOffset(0.9)
        eft_hist.GetXaxis().SetRangeUser(0,2.900)
        eft_hist.GetXaxis().SetLabelSize(0.055)
        eft_hist.GetXaxis().SetTitleSize(0.06)
        eft_hist.GetXaxis().SetTitleOffset(1.2)
        
        if not os.path.isdir(args.InputDir + "/Kinematics_SMvsEFTvsInterference"): os.mkdir(args.InputDir + "/Kinematics_SMvsEFTvsInterference")
    
        c.SaveAs(args.InputDir + "/Kinematics_SMvsEFTvsInterference/SMvsEFTvsInterference_%s.png"%b)
        c.SaveAs(args.InputDir + "/Kinematics_SMvsEFTvsInterference/SMvsEFTvsInterference_%s.pdf"%b)
        

sys.exit(1)
    
colors = [i+5 for i in range(len(dict))]
hist_dict = {}
for name,chain in dict.iteritems():
    print "h_"+name
    hist_dict.update({name:TH1D("h_"+name,"",nbins, xmin,xmax)})   
    chain.Draw( variable_name + " >> h_"+name )
    print hist_dict[name].GetEntries()

# Draw the histograms on a canvas
gStyle.SetOptStat(0) # uncomment if you want to hide statistics box (# events, mean, ...)
c = TCanvas("c","c",600,600)

l = TLegend(0.2,0.7,0.9,0.88)
l.SetBorderSize(0)
l.SetFillStyle(0)
l.SetNColumns(4)
c.cd()
c.SetLogy(0)
c.SetLeftMargin(0.15)
c.SetBottomMargin(0.15)

counter=0
for name,hist in hist_dict.iteritems():
    if hist.GetEntries() != 0: hist.Scale(1./hist.Integral())
    hist.SetLineStyle(colors_dict[name][0])
    hist.SetLineColor(colors_dict[name][1])
    hist.SetLineWidth(colors_dict[name][2])
    hist.GetXaxis().SetTitle(variable_name)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitle("entries (normalized)")
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleOffset(1.1)
    #hist.GetYaxis().SetRangeUser(0.0001,10)
    hist.GetYaxis().SetRangeUser(0.0001,0.4)
    if counter == 0: hist.Draw("hist")
    else: hist.Draw("hist same")
    l.AddEntry(hist,name,"l")
    counter = counter + 1

l.Draw("same")

# Save the Canvas
if not os.path.isdir(output_dir): os.mkdir(output_dir)
c.SaveAs(output_dir + "/" + variable_name + ".pdf")
c.SaveAs(output_dir + "/" + variable_name + ".png")





def main():
    
    gROOT.SetBatch(1)

    parser = ArgumentParser()
    #parser.add_argument('--ncpu', type=int, default=-1,help='number of CPU to use in parallel')
    #parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_Wed15Nov2017_15h36m34s", help='path to the converted delphes files')
    parser.add_argument('--InputDir', default = "CONVERTED_DELPHES_ttbb_ForValidationPerCouplings_fixed", help='path to the converted delphes files')
    parser.add_argument('--coupling', default="cQb8",help='name of coupling')
    args = parser.parse_args()

    files = [i for i in os.listdir(args.InputDir) if ".root" in i if args.coupling in i]
    
    couplings = ["SM","cQQ1","cQQ8","cQt1","cQb1","ctb1","cQt8","cQb8","ctb8"]
    #couplings = ["SM", "C83Qq", "C81Qq", "C13Qq", "C11Qq", "C8td", "C8tu", "C1td", "C1tu", "C8Qd", "C8Qu", "C8tq", "C1Qd", "C1Qu", "C1tq"]
    # for LHCTOPWG validation
    #couplings = ["SM", "cQq83", "cQq81", "cQq13", "cQq11", "ctd8", "ctu8", "ctd1", "ctu1", "cQd8", "cQu8", "ctq8", "cQd1", "cQu1", "ctq1"]
    
    chain_dict = {}
    for c in couplings:
        chain_dict.update({c:TChain("tree")})
    #SM_chain = TChain("")
    #C1tu_chain = TChain("")
    print chain_dict

    for f in files:
        if "SM" in f: chain_dict["SM"].Add(args.InputDir + "/" + f)
        elif "EFT" in f:
            coupling_name = f.split("_")[0]#[:-3]
            #coupling_name = f.split("_")[0]
            chain_dict[coupling_name].Add(args.InputDir + "/" + f)

    branchnames = [i.GetName() for i in chain_dict["SM"].GetListOfBranches()]
    for b in branchnames:
        xmin = min([i.GetMinimum(b) for name,i in chain_dict.iteritems()])
        xmax = max([i.GetMaximum(b) for name,i in chain_dict.iteritems()])
        print b, xmin, xmax
        DrawFromFile(chain_dict, b, 30, xmin,xmax,args.InputDir+"/histograms")


    
if __name__ == "__main__":
    main()