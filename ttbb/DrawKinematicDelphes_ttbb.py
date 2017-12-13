from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TLegend, gROOT, gStyle
import sys
import ROOT
import os
import time
from argparse import ArgumentParser
from array import array
from math import *
import numpy as np



def DrawFromFile(dict, variable_name, nbins, xmin,xmax,output_dir):
   #  colors_dict = { #[color,style,width]
#         'SM': [1,1,3],
#         #LLLL
#         'C11Qq': [1,4,2],
#         'C81Qq': [1,2,2],
#         'C13Qq': [1,8,2],
#         'C83Qq': [1,6,2],
#         #RRRR
#         'C1tu': [3,4,2],
#         'C8tu': [3,2,2],
#         'C1td': [3,8,2],
#         'C8td': [3,6,2],
#          #LLRR
#         'C8tq': [7,2,2],
#         'C8Qu': [7,28,2],
#         'C8Qd': [7,6,2],
#         'C1tq': [7,4,2],
#         'C1Qu': [7,7,2],
#         'C1Qd': [7,8,2]    
#     }
    colors_dict = { #[color,style,width]
        'SM': [1,1,3],
        #LLLL
        "cQQ1": [1,4,2],
        "cQQ8": [1,2,2],
        #LLRR
        "cQt1": [3,6,2],
        "cQb1":[7,6,2],
        "cQt8": [3,8,2],
        "cQb8": [7,8,2],
        #RRRR
        "ctb1": [7,28,2],
        "ctb8": [7,7,2]   
    }
    # For LHCTOPWG validation:
    # colors_dict = { #[color,style,width]
#         'SM': [1,1,2],
#         #LLLL
#         'cQq11': [1,4,1],
#         'cQq81': [1,2,1],
#         'cQq13': [1,8,1],
#         'cQq83': [1,6,1],
#         #RRRR
#         'ctu1': [3,4,1],
#         'ctu8': [3,2,1],
#         'ctd1': [3,8,1],
#         'ctd8': [3,6,1],
#          #LLRR
#         'ctq8': [7,2,1],
#         'cQu8': [7,28,1],
#         'cQd8': [7,6,1],
#         'ctq1': [7,4,1],
#         'cQu1': [7,7,1],
#         'cQd1': [7,8,1]    
#     }
    
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
    c.SetLogy(1)
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
        hist.GetYaxis().SetRangeUser(0.0001,10)
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
    parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_TRAINING", help='path to the converted delphes files')
    #parser.add_argument('--tag', default=time.strftime("%a%d%b%Y_%Hh%Mm%Ss"),help='name of output directory')
    args = parser.parse_args()

    files = [i for i in os.listdir(args.InputDir) if ".root" in i]
    
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