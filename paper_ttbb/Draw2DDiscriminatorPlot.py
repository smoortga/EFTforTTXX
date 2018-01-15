import os
from argparse import ArgumentParser
import sys
import matplotlib.pyplot as plt
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

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_ForValidationPerCouplings_fixed", help='path to the directory of all restricted processes')
parser.add_argument('--TrainingFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/model_checkpoint_save.hdf5", help='path to the directory of SM-like xsec measurement')
parser.add_argument('--ScalerFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Training_fixed/training_output/scaler.pkl", help='path to the directory of SM-like xsec measurement')
parser.add_argument('--pickEvery', type=int, default=1, help='pick one event every so many events')

args = parser.parse_args()


def GetDiscrSMvsEFT(filename,model,scaler):
    X = rootnp.root2array(args.ValidationDir + "/" + filename,"tree")
    X = rootnp.rec2array(X)
    X = scaler.transform(X)
    discr = model.predict(X)
    return [i+j for i,j in zip(discr[:,1],discr[:,2])]

def GetDiscrtLvstR(filename,model,scaler):
    X = rootnp.root2array(args.ValidationDir + "/" + filename,"tree")
    X = rootnp.rec2array(X)
    X = scaler.transform(X)
    discr = model.predict(X)
    return [i/(i+j) for i,j in zip(discr[:,1],discr[:,2])]
    



files_dict = {
    "SM":["cQb8_0p0_run_01_tag_1_delphes_events.root","cQb1_0p0_run_01_tag_1_delphes_events.root","cQQ1_0p0_run_01_tag_1_delphes_events.root","cQQ8_0p0_run_01_tag_1_delphes_events.root"],
    "cQb1_20p0":["cQb1_20p0_run_01_tag_1_delphes_events.root","cQb8_20p0_run_01_tag_1_delphes_events.root","cQQ8_20p0_run_01_tag_1_delphes_events.root","cQQ1_20p0_run_01_tag_1_delphes_events.root"],
    "ctb1_20p0":["ctb1_20p0_run_01_tag_1_delphes_events.root","ctb1_20p0_run_01_tag_1_delphes_events.root","cQt1_20p0_run_01_tag_1_delphes_events.root","cQt8_20p0_run_01_tag_1_delphes_events.root"]
}

SMvsEFT_discr_dict = {}
tLvstR_discr_dict = {}
model = load_model(args.TrainingFile)
scaler = pickle.load(open(args.ScalerFile,'r'))
for c,files in files_dict.iteritems():
    for idx,name in enumerate(files):
        print "Processing: %s"%name
        if idx == 0:
            SMvsEFT_discr_dict[c] = GetDiscrSMvsEFT(name,model,scaler)
            tLvstR_discr_dict[c] = GetDiscrtLvstR(name,model,scaler)
        else:
            SMvsEFT_discr_dict[c] = np.concatenate((SMvsEFT_discr_dict[c],GetDiscrSMvsEFT(name,model,scaler)))
            tLvstR_discr_dict[c] = np.concatenate((tLvstR_discr_dict[c],GetDiscrtLvstR(name,model,scaler)))

max_len = min(len(SMvsEFT_discr_dict["SM"]),len(SMvsEFT_discr_dict["cQb1_20p0"]),len(SMvsEFT_discr_dict["ctb1_20p0"]))
SMvsEFT_discr_dict["SM"] = SMvsEFT_discr_dict["SM"][0:max_len]
SMvsEFT_discr_dict["cQb1_20p0"] = SMvsEFT_discr_dict["cQb1_20p0"][0:max_len]
SMvsEFT_discr_dict["ctb1_20p0"] = SMvsEFT_discr_dict["ctb1_20p0"][0:max_len]
tLvstR_discr_dict["SM"] = tLvstR_discr_dict["SM"][0:max_len]
tLvstR_discr_dict["cQb1_20p0"] = tLvstR_discr_dict["cQb1_20p0"][0:max_len]
tLvstR_discr_dict["ctb1_20p0"] = tLvstR_discr_dict["ctb1_20p0"][0:max_len]

ROOT.gROOT.SetBatch(True)

# PREPARE ALL FILES


# Eff_vs_Disc_2D_SM = ROOT.TGraph2D()
# Eff_vs_Disc_2D_EFT_tL = ROOT.TGraph2D()
# Eff_vs_Disc_2D_EFT_tR = ROOT.TGraph2D()

disc_min = 0
disc_max = 1
ncuts = 25

# SMvsEFT_Disc = np.arange(disc_min,disc_max,(disc_max-disc_min)/float(ncuts))
# tLvstR_Disc = np.arange(disc_min,disc_max,(disc_max-disc_min)/float(ncuts))

Plot_2D_SM = ROOT.TH2D("Plot_2D_SM"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/2.)
Plot_2D_EFT_tL = ROOT.TH2D("Plot_2D_EFT_tL"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/2.)
Plot_2D_EFT_tR = ROOT.TH2D("Plot_2D_EFT_tR"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/2.)


# SMvsEFT_Disc_bin = []
# tLvstR_Disc_bin = []
# for idx,i in enumerate(SMvsEFT_Disc):
# 	SMvsEFT_Disc_bin.insert(idx,Plot_2D_SM.GetXaxis().FindBin(i))
# 	tLvstR_Disc_bin.insert(idx,Plot_2D_SM.GetXaxis().FindBin(i))
# 	
	
# FOR SM

for i in xrange(len(SMvsEFT_discr_dict["SM"])):
	if (i%args.pickEvery)!=0:continue
	Plot_2D_SM.Fill(SMvsEFT_discr_dict["SM"][i],tLvstR_discr_dict["SM"][i])

# Total_Integral_SM = Plot_2D_SM.Integral(0,ncuts,0,ncuts)
# for idx,i in enumerate(SMvsEFT_Disc):
# 	for jdx,j in enumerate(tLvstR_Disc):
# 		Eff_vs_Disc_2D_SM.SetPoint(jdx+(ncuts*idx),i,j,Plot_2D_SM.Integral(SMvsEFT_Disc_bin[idx],ncuts,tLvstR_Disc_bin[jdx],ncuts) / Total_Integral_SM)

print "Done filling histograms for SM"


# FOR EFT tL

for i in xrange(len(SMvsEFT_discr_dict["cQb1_20p0"])):
	if (i%args.pickEvery)!=0:continue
	Plot_2D_EFT_tL.Fill(SMvsEFT_discr_dict["cQb1_20p0"][i],tLvstR_discr_dict["cQb1_20p0"][i])

# Total_Integral_EFT_tL = Plot_2D_EFT_tL.Integral(0,ncuts,0,ncuts)
# for idx,i in enumerate(SMvsEFT_Disc):
# 	for jdx,j in enumerate(tLvstR_Disc):
# 		Eff_vs_Disc_2D_EFT_tL.SetPoint(jdx+(ncuts*idx),i,j,Plot_2D_EFT_tL.Integral(SMvsEFT_Disc_bin[idx],ncuts,tLvstR_Disc_bin[jdx],ncuts) / Total_Integral_EFT_tL)

print "Done filling histograms for EFT tL"


# FOR EFT tR

for i in xrange(len(SMvsEFT_discr_dict["ctb1_20p0"])):
	if (i%args.pickEvery)!=0:continue
	Plot_2D_EFT_tR.Fill(SMvsEFT_discr_dict["ctb1_20p0"][i],tLvstR_discr_dict["ctb1_20p0"][i])

# Total_Integral_EFT_tR = Plot_2D_EFT_tR.Integral(0,ncuts,0,ncuts)
# for idx,i in enumerate(SMvsEFT_Disc):
# 	for jdx,j in enumerate(tLvstR_Disc):
# 		Eff_vs_Disc_2D_EFT_tR.SetPoint(jdx+(ncuts*idx),i,j,Plot_2D_EFT_tR.Integral(SMvsEFT_Disc_bin[idx],ncuts,tLvstR_Disc_bin[jdx],ncuts) / Total_Integral_EFT_tR)


print "Done filling histograms for EFT tR"






Plot_2D_SM.Scale(1/Plot_2D_SM.Integral())
Plot_2D_EFT_tL.Scale(1/Plot_2D_EFT_tL.Integral())
Plot_2D_EFT_tR.Scale(1/Plot_2D_EFT_tR.Integral())


Plot_2D_SM_most = ROOT.TH2D("Plot_2D_SM_most"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_SM_middle = ROOT.TH2D("Plot_2D_SM_middle"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_SM_least = ROOT.TH2D("Plot_2D_SM_least"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tL_most = ROOT.TH2D("Plot_2D_EFT_tL_most"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tL_middle = ROOT.TH2D("Plot_2D_EFT_tL_middle"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tL_least = ROOT.TH2D("Plot_2D_EFT_tL_least"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tR_most = ROOT.TH2D("Plot_2D_EFT_tR_most"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tR_middle = ROOT.TH2D("Plot_2D_EFT_tR_middle"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)
Plot_2D_EFT_tR_least = ROOT.TH2D("Plot_2D_EFT_tR_least"," ;discriminator P(t_{L}) + P(t_{R});discriminator #frac{P(t_{L})}{P(t_{L}) + P(t_{R})}",ncuts,disc_min,disc_max,ncuts,disc_min,disc_max+(disc_max-disc_min)/10.)

for i in range(Plot_2D_SM.GetNbinsX()+1):
    for j in range(Plot_2D_SM.GetNbinsY()):
        B_content =Plot_2D_SM.GetBinContent(i,j)
        C_content =Plot_2D_EFT_tL.GetBinContent(i,j)
        DUSG_content =Plot_2D_EFT_tR.GetBinContent(i,j)
        if B_content == 0 and C_content == 0 and DUSG_content == 0: continue
        elif B_content > C_content and B_content > DUSG_content:
            Plot_2D_SM_most.SetBinContent(i,j,B_content)
            if C_content > DUSG_content: 
                Plot_2D_EFT_tL_middle.SetBinContent(i,j,C_content)
                Plot_2D_EFT_tR_least.SetBinContent(i,j,DUSG_content)
            else: 
                Plot_2D_EFT_tR_middle.SetBinContent(i,j,DUSG_content)
                Plot_2D_EFT_tL_least.SetBinContent(i,j,C_content)
        elif C_content > B_content and C_content > DUSG_content:
            Plot_2D_EFT_tL_most.SetBinContent(i,j,C_content)
            if B_content > DUSG_content: 
                Plot_2D_SM_middle.SetBinContent(i,j,B_content)
                Plot_2D_EFT_tR_least.SetBinContent(i,j,DUSG_content)
            else: 
                Plot_2D_EFT_tR_middle.SetBinContent(i,j,DUSG_content)
                Plot_2D_SM_least.SetBinContent(i,j,B_content)
        elif DUSG_content > C_content and DUSG_content > B_content:
            Plot_2D_EFT_tR_most.SetBinContent(i,j,DUSG_content)
            if C_content > B_content: 
                Plot_2D_EFT_tL_middle.SetBinContent(i,j,C_content)
                Plot_2D_SM_least.SetBinContent(i,j,B_content)
            else: 
                Plot_2D_SM_middle.SetBinContent(i,j,B_content)
                Plot_2D_EFT_tL_least.SetBinContent(i,j,C_content)

Plot_2D_SM_most.SetFillColor(ROOT.kRed)
Plot_2D_SM_middle.SetFillColorAlpha(ROOT.kRed,0.6)
Plot_2D_SM_least.SetFillColorAlpha(ROOT.kRed,0.3)
Plot_2D_EFT_tR_most.SetFillColor(ROOT.kBlue)
Plot_2D_EFT_tR_middle.SetFillColorAlpha(ROOT.kBlue,0.6)
Plot_2D_EFT_tR_least.SetFillColorAlpha(ROOT.kBlue,0.3)
Plot_2D_EFT_tL_most.SetFillColor(8)
Plot_2D_EFT_tL_middle.SetFillColorAlpha(8,0.8)
Plot_2D_EFT_tL_least.SetFillColorAlpha(8,0.5)

# Plot_2D_SM_middle.SetFillStyle(3002)
# Plot_2D_SM_least.SetFillStyle(3003)
# Plot_2D_EFT_tR_middle.SetFillStyle(3002)
# Plot_2D_EFT_tR_least.SetFillStyle(3003)
# Plot_2D_EFT_tL_middle.SetFillStyle(3002)
# Plot_2D_EFT_tL_least.SetFillStyle(3003)

hs = ROOT.THStack("hs","")
hs.Add(Plot_2D_SM_most)
hs.Add(Plot_2D_EFT_tL_most)
hs.Add(Plot_2D_EFT_tR_most)
hs.Add(Plot_2D_SM_middle)
hs.Add(Plot_2D_EFT_tL_middle)
hs.Add(Plot_2D_EFT_tR_middle)
hs.Add(Plot_2D_SM_least)
hs.Add(Plot_2D_EFT_tL_least)
hs.Add(Plot_2D_EFT_tR_least)

min_z = min([Plot_2D_SM_most.GetMinimum(),Plot_2D_EFT_tL_most.GetMinimum(),Plot_2D_EFT_tR_most.GetMinimum()])
max_z = max([Plot_2D_SM_most.GetMaximum(),Plot_2D_EFT_tL_most.GetMaximum(),Plot_2D_EFT_tR_most.GetMaximum()])
Plot_2D_SM_most.SetMinimum(min_z)
Plot_2D_SM_most.SetMaximum(max_z)
Plot_2D_SM_middle.SetMinimum(min_z)
Plot_2D_SM_middle.SetMaximum(max_z)
Plot_2D_SM_least.SetMinimum(min_z)
Plot_2D_SM_least.SetMaximum(max_z)
Plot_2D_EFT_tL_most.SetMinimum(min_z)
Plot_2D_EFT_tL_most.SetMaximum(max_z)
Plot_2D_EFT_tL_middle.SetMinimum(min_z)
Plot_2D_EFT_tL_middle.SetMaximum(max_z)
Plot_2D_EFT_tL_least.SetMinimum(min_z)
Plot_2D_EFT_tL_least.SetMaximum(max_z)
Plot_2D_EFT_tR_most.SetMinimum(min_z)
Plot_2D_EFT_tR_most.SetMaximum(max_z)
Plot_2D_EFT_tR_middle.SetMinimum(min_z)
Plot_2D_EFT_tR_middle.SetMaximum(max_z)
Plot_2D_EFT_tR_least.SetMinimum(min_z)
Plot_2D_EFT_tR_least.SetMaximum(max_z)

## DRAW 2D DISCRIMINATOR DISTRIBUTIONS AND SCATTER PLOT

c_sc = ROOT.TCanvas("c_sc","c_sc",800,600)
c_sc.cd()
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetMargin(0.2,0.1,0.15,0.07)
l = ROOT.TLegend(0.23,0.75,0.90,0.85)
l.SetNColumns(3)
l.SetFillColor(0)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.SetTextSize(0.045)


Plot_2D_SM_most.GetYaxis().SetTitleOffset(1.3)
Plot_2D_SM_most.GetYaxis().SetTitleSize(0.055)
Plot_2D_SM_most.GetYaxis().SetLabelSize(0.045)
#Plot_2D_SM_most.GetYaxis().SetRangeUser(0,1.5)
Plot_2D_SM_most.GetXaxis().SetTitleOffset(1.2)
Plot_2D_SM_most.GetXaxis().SetTitleSize(0.055)
Plot_2D_SM_most.GetXaxis().SetLabelSize(0.045)
Plot_2D_SM_most.GetZaxis().SetLabelSize(0.055)
Plot_2D_EFT_tL.GetYaxis().SetTitleOffset(1.1)
Plot_2D_EFT_tL.GetYaxis().SetTitleSize(0.055)
Plot_2D_EFT_tL.GetYaxis().SetLabelSize(0.045)
#Plot_2D_EFT_tL.GetYaxis().SetRangeUser(0,1.5)
Plot_2D_EFT_tL.GetXaxis().SetTitleSize(0.055)
Plot_2D_EFT_tL.GetXaxis().SetLabelSize(0.055)
Plot_2D_EFT_tL.GetZaxis().SetLabelSize(0.055)
Plot_2D_EFT_tR.GetYaxis().SetTitleOffset(1.1)
Plot_2D_EFT_tR.GetYaxis().SetTitleSize(0.055)
Plot_2D_EFT_tR.GetYaxis().SetLabelSize(0.045)
#Plot_2D_EFT_tR.GetYaxis().SetRangeUser(0,1.5)
Plot_2D_EFT_tR.GetXaxis().SetTitleSize(0.065)
Plot_2D_EFT_tR.GetXaxis().SetLabelSize(0.055)
Plot_2D_EFT_tR.GetZaxis().SetLabelSize(0.055)

Plot_2D_SM.SetFillColor(ROOT.kRed)
Plot_2D_EFT_tL.SetFillColor(8)
Plot_2D_EFT_tR.SetFillColor(ROOT.kBlue)
l.AddEntry(Plot_2D_SM,"SM","f")
l.AddEntry(Plot_2D_EFT_tL,"#splitline{SM+EFT}{left-handed top}","f")
l.AddEntry(Plot_2D_EFT_tR,"#splitline{SM+EFT}{right-handed top}","f")

c_sc.cd()
ROOT.gPad.SetTickx(1)
ROOT.gPad.SetTicky(1)	
Plot_2D_SM_most.Draw("box")
Plot_2D_EFT_tL_most.Draw("box same")
Plot_2D_EFT_tR_most.Draw("box same")
Plot_2D_SM_middle.Draw("box same")
Plot_2D_EFT_tL_middle.Draw("box same")
Plot_2D_EFT_tR_middle.Draw("box same")
Plot_2D_SM_least.Draw("box same")
Plot_2D_EFT_tL_least.Draw("box same")
Plot_2D_EFT_tR_least.Draw("box same")

l.Draw("same")

if not os.path.isdir(args.InputDir + "/2D_Discriminator_Plots"): os.mkdir(args.InputDir + "/2D_Discriminator_Plots")
c_sc.SaveAs(args.InputDir + "/2D_Discriminator_Plots/2D_Plot_Box.png")
c_sc.SaveAs(args.InputDir + "/2D_Discriminator_Plots/2D_Plot_Box.pdf")