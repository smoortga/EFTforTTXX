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
    colors_dict = { #[color,style,width]
        'SM': [1,1,2],
        #LLLL
        'C11Qq': [1,4,1],
        'C81Qq': [1,2,1],
        'C13Qq': [1,8,1],
        'C83Qq': [1,6,1],
        #RRRR
        'C1tu': [3,4,1],
        'C8tu': [3,2,1],
        'C1td': [3,8,1],
        'C8td': [3,6,1],
         #LLRR
        'C8tq': [7,2,1],
        'C8Qu': [7,28,1],
        'C8Qd': [7,6,1],
        'C1tq': [7,4,1],
        'C1Qu': [7,7,1],
        'C1Qd': [7,8,1]    
    }
    
    colors = [i+5 for i in range(len(dict))]
    hist_dict = {}
    for name,chain in dict.iteritems():
        hist_dict.update({name:TH1D("h_"+name,"",nbins, xmin,xmax)})   
        chain.Draw( variable_name + " >> h_"+name )

    # Draw the histograms on a canvas
    gStyle.SetOptStat(0) # uncomment if you want to hide statistics box (# events, mean, ...)
    c = TCanvas("c","c",600,500)

    l = TLegend(0.6,0.5,0.9,0.88)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetNColumns(2)
    c.cd()
    c.SetLogy(1)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)

    counter=0
    for name,hist in hist_dict.iteritems():
        hist.Scale(1./hist.Integral())
        hist.SetLineStyle(colors_dict[name][0])
        hist.SetLineColor(colors_dict[name][1])
        hist.SetLineWidth(colors_dict[name][2])
        hist.GetXaxis().SetTitle(variable_name)
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitle("entries (normalized)")
        hist.GetYaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetRangeUser(0.00001,1)
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
    parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_DILEPTON", help='path to the converted delphes files')
    #parser.add_argument('--tag', default=time.strftime("%a%d%b%Y_%Hh%Mm%Ss"),help='name of output directory')
    args = parser.parse_args()

    files = [i for i in os.listdir(args.InputDir) if ".root" in i]
    
    couplings = ["SM", "C83Qq", "C81Qq", "C13Qq", "C11Qq", "C8td", "C8tu", "C1td", "C1tu", "C8Qd", "C8Qu", "C8tq", "C1Qd", "C1Qu", "C1tq"]
    
    chain_dict = {}
    for c in couplings:
        chain_dict.update({c:TChain("tree")})
    #SM_chain = TChain("")
    #C1tu_chain = TChain("")
    print chain_dict

    for f in files:
        if "Order0" in f: chain_dict["SM"].Add(args.InputDir + "/" + f)
        elif "Order2" in f:
            coupling_name = f.split("_")[0][:-3]
            chain_dict[coupling_name].Add(args.InputDir + "/" + f)

    branchnames = [i.GetName() for i in chain_dict["SM"].GetListOfBranches()]
    for b in branchnames:
        xmin = min([i.GetMinimum(b) for name,i in chain_dict.iteritems()])
        xmax = max([i.GetMaximum(b) for name,i in chain_dict.iteritems()])
        print b
        DrawFromFile(chain_dict, b, 30, xmin,xmax,args.InputDir+"/histograms")


    
if __name__ == "__main__":
    main()