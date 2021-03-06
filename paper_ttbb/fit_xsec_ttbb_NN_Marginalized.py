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
import math
from array import array
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans-serif')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttbb_DiLepton_ForValidationPerCouplings_Marginalized", help='path to the directory of all restricted processes')
parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_ttbb_Marginalized", help='path to the directory of all restricted processes')
#parser.add_argument('--SMxsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0_backup/MODEL_ttcc_inclusive_SMRestrictedMassAndCabibbo", help='path to the directory of SM-like xsec measurement')
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

def NN_validate(filename,cutEFT, cuttL, cuttR,original_n_events = 30000):
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
    discr_dict = {}
    for class_n in set(i for j,i in classes_dict.iteritems()):
        discr_dict[class_n] = model.predict(X)[:,class_n]
    discr_SMvsEFT = np.asarray([(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1])])
    discr_tLvstR = np.asarray([discr_dict[1][jdx]/(discr_dict[1][jdx]+discr_dict[2][jdx]) for jdx,j in enumerate(discr_dict[1])])
    nEvents = len(discr_SMvsEFT)
    print float(len(discr_SMvsEFT)), sum_original_n_events, 100*float(len(discr_SMvsEFT))/float(original_n_events),"%"

    zipped = zip(discr_SMvsEFT,discr_tLvstR)
    discr_EFT = [i for i in zipped if i[0]>cutEFT[0] and i[1] > cutEFT[1]]
    discr_tL = [i for i in zipped if i[0]>cuttL[0] and i[1] > cuttL[1]]
    discr_tR = [i for i in zipped if i[0]>cuttR[0] and i[1] < cuttR[1]]
    
    print "selection efficiency EFT NN cut: " ,100*float(len(discr_EFT))/float(nEvents)
    print "selection efficiency tL NN cut: " ,100*float(len(discr_tL))/float(nEvents)
    print "selection efficiency tR NN cut: " ,100*float(len(discr_tR))/float(nEvents)

    #print ""
    #print filename, float(len(discr)),"/",float(len(filename)*original_n_events),"%"
    return float(len(discr_EFT))/float(original_n_events),float(len(discr_tL))/float(original_n_events),float(len(discr_tR))/float(original_n_events)


def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    # start with first coupling
    c1 = coupling.split("_")[0]
    val1 = val.split("_")[0] + "_" + val.split("_")[1]
    res1=-999
    if "minus" in val1: res1= -float(val1.replace(c1+"_minus","").replace("p","."))
    else: res1= float(val1.replace(c1+"_","").replace("p","."))
    # now the second
    c2 = coupling.split("_")[1]
    val2 = val.split("_")[2] + "_" + val.split("_")[3]
    res2=-999
    if "minus" in val2: res2= -float(val2.replace(c2+"_minus","").replace("p","."))
    else: res2= float(val2.replace(c2+"_","").replace("p","."))
    return [res1,res2]
    
def extract_coupling_string(coupling,val):
    """ extract the string value of a coupling in the directory name"""
    return val.split("_")[1] + "_" + val.split("_")[3]
    
def func2(x, a, b, c):
        return (a + b*x + c*x*x)   

def func4(x, a, b, c, d, e):
        return (a + b*x + c*x*x + d*x*x*x + e*x*x*x*x)    

def func(x, p0, p1a, p1b, p2aa, p2bb, p2ab):
        return p0 + p1a*x[0] + p1b*x[1] + p2aa*x[0]*x[0] + p2bb*x[1]*x[1] + p2ab*x[0]*x[1]
    

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d and not 'Limits' in d and not "Discriminator" in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]

validation_files = {}
files = [f for f in os.listdir(args.ValidationDir) if ".root" in f]
couplings = set([f.split("_")[0] + "_" + f.split("_")[2]  for f in os.listdir(args.ValidationDir) if ".root" in f])
coupl_strengths = set([f.split("_")[1] + "_" + f.split("_")[3] for f in os.listdir(args.ValidationDir) if ".root" in f])
for c in couplings:
    validation_files[c]={}
    for s in coupl_strengths:
        validation_files[c][s]=[]
for f in files:
    coupling_name = f.split("_")[0] + "_" + f.split("_")[2]
    strength_name = f.split("_")[1] + "_" + f.split("_")[3]
    validation_files[coupling_name][strength_name].append(f)


# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1a = []
p1ae = []
p1b = []
p1be = []
p2aa = []
p2aae = []
p2bb = []
p2bbe = []
p2ab = []
p2abe = []
coupling_names = []


nevents={"run_01":30000,"run_02":0}

result_dict = {}

EFT_cuts = [0.45,0]
tR_cuts = [0.45,0.5]
tL_cuts = [0.45,0.5]

# for coupling,orders in process_dict.iteritems():
#     result_dict[coupling] = {}
#     for o in sorted(orders):
#         #if extract_coupling(coupling,o)[0] != 0.0 or extract_coupling(coupling,o)[1] != 0.0: continue
#         
#         print "Processing %s %.1f,%.1f"%(coupling,extract_coupling(coupling,o)[0],extract_coupling(coupling,o)[1])
#         ##################
#         #
#         # Get cross section
#         #
#         ##################
#         if os.path.isfile(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt"):
#             f_ = open(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt", 'r')
#             lines = f_.readlines()
#             xsec = float(lines[-1].split(" ")[2])
#             xsec_error = float(lines[-1].split(" ")[3])
#         else: continue
#         #xsec = 0.01
#         #xsec_error = 0.001
#         
#         #print xsec, xsec_error
#         
#         ##################
#         #
#         # NN output
#         #
#         ##################
#         
#         filelist = validation_files[coupling][extract_coupling_string(coupling,o)]
#         #filelist = validation_files[coupling][o.split("_")[1]]
#         if len(filelist) == 0: continue
#         sum_original_n_events = 0
#         for f in filelist:
#             sum_original_n_events += get_Original_nevents(f,nevents)
#         frac_passing_evts_EFT, frac_passing_evts_tL, frac_passing_evts_tR = NN_validate(filelist,cutEFT=EFT_cuts, cuttL=tL_cuts, cuttR = tR_cuts,original_n_events = sum_original_n_events)
#         result_dict[coupling][o] = [xsec,xsec_error,frac_passing_evts_EFT, frac_passing_evts_tL, frac_passing_evts_tR, sum_original_n_events]
#    
# if not os.path.isdir(args.ValidationDir + "/validation_output"): os.mkdir(args.ValidationDir + "/validation_output")
# pickle.dump(result_dict,open(args.ValidationDir + "/validation_output/result_dict.pkl",'wb'))
result_dict = pickle.load(open(args.ValidationDir + "/validation_output/result_dict.pkl",'r'))

print "done processing NN outputs"


# find out the mean of the SM contributions as an average value
# sm_values = []
# for coupling,orders in result_dict.iteritems():
#     for o,results in orders.iteritems():
#         if o == "0p0": sm_values.append(results[0]*results[2])
# sm_xsec = np.mean(np.asarray(sm_values))
# sm_xsec_error = np.std(np.asarray(sm_values))
    
# sm_xsec = 0.085*0.11*0.05 # pb from TOP-16-010
# sm_xsec_error = 0.041*0.11*0.05# pb from TOP-16-010
sm_xsec_frac_error = 0.1
sm_xsec_frac_error_300fb = 0.1

fracerr_measured_sm = sm_xsec_frac_error

# define a file for the limits to be stored: 
# limits_file = open('%s/limits_ttbb_cutonNN.txt'%(args.InputDir),"w")
# #limits_file = open('%s/limits_ttbb_NN_Binary_cut0p75.txt'%(args.InputDir),"w")
# #limits_file = open('%s/limits_ttbb_nocut.txt'%(args.InputDir),"w")
# limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))


for coupling,values in result_dict.iteritems():
    if "Discriminator" in coupling: continue
    #if coupling == "cQQ1": continue
    coupling_strengths_x1 = []
    coupling_strengths_x2 = []
    xsec_EFT = []
    xsec_EFT_error = []
    xsec_tL = []
    xsec_tL_error = []
    xsec_tR = []
    xsec_tR_error = []
    #nevents = []
    sm_xsec_EFT = 0
    sm_xsec_EFT_error = 0
    sm_xsec_tL = 0
    sm_xsec_tL_error = 0
    sm_xsec_tR = 0
    sm_xsec_tR_error = 0

    for v,result in values.iteritems():
        #print v
        if v == "cQb1_0p0_ctb1_0p0": 
            sm_xsec_EFT = result[0]*result[2]
            sm_xsec_EFT_error = math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[5]*result[2])/result[5])**2 )
            sm_xsec_tL = result[0]*result[3]
            sm_xsec_tL_error = math.sqrt( (result[1]*result[3])**2 + (result[0]*math.sqrt(result[5]*result[3])/result[5])**2 )
            sm_xsec_tR = result[0]*result[4]
            sm_xsec_tR_error = math.sqrt( (result[1]*result[4])**2 + (result[0]*math.sqrt(result[5]*result[4])/result[5])**2 )
        coupling_strengths_x1.append(extract_coupling(coupling,v)[0])
        coupling_strengths_x2.append(extract_coupling(coupling,v)[1])
        xsec_EFT.append(result[0]*result[2])
        xsec_EFT_error.append(math.sqrt( (result[1]*result[2])**2 + (result[0]*math.sqrt(result[5]*result[2])/result[5])**2 ))
        xsec_tL.append(result[0]*result[3])
        xsec_tL_error.append(math.sqrt( (result[1]*result[3])**2 + (result[0]*math.sqrt(result[5]*result[3])/result[5])**2 ))
        xsec_tR.append(result[0]*result[4])
        xsec_tR_error.append(math.sqrt( (result[1]*result[4])**2 + (result[0]*math.sqrt(result[5]*result[4])/result[5])**2 ))


    #print xsec,xsec_error
    coeff_EFT, covmat_EFT = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),xsec_EFT,sigma=xsec_EFT_error)#absolute_sigma
    errors_EFT = sqrt(diag(covmat_EFT))
    
    coeff_tL, covmat_tL = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),xsec_tL,sigma=xsec_tL_error)#absolute_sigma
    errors_tL = sqrt(diag(covmat_tL))
    
    coeff_tR, covmat_tR = curve_fit(func,np.vstack((coupling_strengths_x1, coupling_strengths_x2)),xsec_tR,sigma=xsec_tR_error)#absolute_sigma
    errors_tR = sqrt(diag(covmat_tR))

    
    
    #**************************
    #
    #   PLOTTING
    #
    #**************************

    x1min = -7
    x1max = +7
    x2min = -7
    x2max = +7
    nsteps = 100
    x1 = np.arange(x1min, x1max+ float(x1max-x1min)/float(nsteps), float(x1max-x1min)/float(nsteps))
    x2 = np.arange(x2min, x2max+ float(x2max-x2min)/float(nsteps), float(x2max-x2min)/float(nsteps))
    side_x1 = np.linspace(x1min, x1max+ float(x1max-x1min)/float(nsteps), nsteps)
    side_x2 = np.linspace(x2min, x2max+ float(x2max-x2min)/float(nsteps), nsteps)
    X1, X2 = np.meshgrid(side_x1, side_x2)
    size = X1.shape
    x1_1d = X1.reshape((1, np.prod(size)))
    x2_1d = X2.reshape((1, np.prod(size)))
    xdata = np.vstack((x1_1d, x2_1d))
    
    y_EFT_down = func(xdata,*(coeff_EFT - [i for i in errors_EFT]))
    y_EFT_tmp = np.asarray([(i-sm_xsec_EFT)**2/((sm_xsec_frac_error*sm_xsec_EFT)**2) for i in y_EFT_down])
    Y_EFT = y_EFT_tmp.reshape(size)
    
    y_tL_down = func(xdata,*(coeff_tL - [i for i in errors_tL]))
    y_tL_tmp = np.asarray([(i-sm_xsec_tL)**2/((sm_xsec_frac_error*sm_xsec_tL)**2) for i in y_tL_down])
    Y_tL = y_tL_tmp.reshape(size)
    
    y_tR_down = func(xdata,*(coeff_tR - [i for i in errors_tR]))
    y_tR_tmp = np.asarray([(i-sm_xsec_tR)**2/((sm_xsec_frac_error*sm_xsec_tR)**2) for i in y_tR_down])
    Y_tR = y_tR_tmp.reshape(size)
    #print sm_xsec_EFT,sm_xsec_frac_error,X1, X2, Y
    # y_up = func(x,*(coeff + [i for i in errors]))
#     y_down = func(x,*(coeff - [i for i in errors]))
    
    y_tmp = np.asarray([i+j for i,j in zip(y_tL_tmp,y_tR_tmp)])
    Y = y_tmp.reshape(size)
    
    
    
    
    levels = [3.84]
    levels_2dof = [5.991]
    CS_EFT = plt.contour(X1, X2, Y_EFT,levels,colors=('r'))
    p_EFT = CS_EFT.collections[0].get_paths()[0]
    v_EFT = p_EFT.vertices
    px_EFT = v_EFT[:,0]
    py_EFT = v_EFT[:,1]
    #plt.clabel(CS_EFT, inline=1, fontsize=10)
    CS_tL = plt.contour(X1, X2, Y_tL,levels,colors=('g'))
    p_tL = CS_tL.collections[0].get_paths()[0]
    v_tL = p_tL.vertices
    px_tL = v_tL[:,0]
    py_tL = v_tL[:,1]
    #plt.clabel(CS_tL, inline=1, fontsize=10)
    CS_tR = plt.contour(X1, X2, Y_tR,levels,colors=('b'))
    p_tR = CS_tR.collections[0].get_paths()[0]
    v_tR = p_tR.vertices
    px_tR = v_tR[:,0]
    py_tR = v_tR[:,1]
    ######
    CS = plt.contour(X1, X2, Y,levels_2dof,colors=('r'),linestyles = 'dashed')
    p = CS.collections[0].get_paths()[0]
    v = p.vertices
    px = v[:,0]
    py = v[:,1]
    


    gr_EFT = ROOT.TGraph(len(px_EFT),array('d',px_EFT),array('d',py_EFT))
    gr_EFT.SetMarkerSize(0)
    gr_EFT.SetFillColor(2)
    gr_EFT.SetLineColor(2)
    gr_EFT.SetLineWidth(2)
    
    gr_tL = ROOT.TGraph(len(px_tL),array('d',px_tL),array('d',py_tL))
    gr_tL.SetMarkerSize(0)
    gr_tL.SetFillColor(8)
    gr_tL.SetLineColor(8)
    gr_tL.SetLineWidth(2)
    
    gr_tR = ROOT.TGraph(len(px_tR),array('d',px_tR),array('d',py_tR))
    gr_tR.SetMarkerSize(0)
    gr_tR.SetFillColor(4)
    gr_tR.SetLineColor(4)
    gr_tR.SetLineWidth(2)
    
    gr_comb = ROOT.TGraph(len(px),array('d',px),array('d',py))
    gr_comb.SetMarkerSize(0)
    gr_comb.SetFillColor(2)
    gr_comb.SetLineColor(2)
    gr_comb.SetLineStyle(2)
    gr_comb.SetLineWidth(2)
    
    mg = ROOT.TMultiGraph("mg",";%s [TeV^{-2}];%s [TeV^{-2}]"%(convertToLatex(coupling.split("_")[0]),convertToLatex(coupling.split("_")[1])))
    mg.Add(gr_EFT,"lc")
    mg.Add(gr_tL,"lc")
    mg.Add(gr_tR,"lc")
    #mg.Add(gr_comb,"lc")
    
    ROOT.TGaxis.SetMaxDigits(2)
    c = ROOT.TCanvas("c","c",800,900)
    ROOT.gPad.SetMargin(0.15,0.05,0.15,0.2)
    #gr_data.Draw("AP")
    mg.Draw("AC")
    #ymin = sm_xsec - 5*sm_xsec_frac_error*sm_xsec
    #ymax = sm_xsec + 10*sm_xsec_frac_error*sm_xsec
    mg.SetMinimum(x2min)
    mg.SetMaximum(x2max)
    mg.GetYaxis().CenterTitle()
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetYaxis().SetTitleOffset(1.05)
    mg.GetYaxis().SetLabelSize(0.05)
    mg.GetXaxis().CenterTitle()
    mg.GetXaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetTitleOffset(1.1)
    mg.GetXaxis().SetLabelSize(0.05)
    mg.GetXaxis().SetLimits(x1min,x1max)
    
    
    l = ROOT.TLegend(0.15,0.8,0.95,0.93)
    l.SetBorderSize(1)
    l.SetNColumns(3)
    l.SetTextSize(0.04)
    l.AddEntry(gr_EFT,"#splitline{SM vs EFT}{95% CL}","l")
    l.AddEntry(gr_tL,"#splitline{SM vs t_{L}}{95% CL}","l")
    #l.AddEntry(gr_comb,"#splitline{combined}{95% CL}","l")
    l.AddEntry(gr_tR,"#splitline{SM vs t_{R}}{95% CL}","l")
    
    
    l.Draw("same")
    
    # text = ROOT.TLatex(0.7,0.9,"M_{4b} < 600 GeV")
#     text.SetTextSize(0.06)
#     text.SetTextFont(42)
#     text.DrawLatexNDC(0.7,0.9,"M_{4b} < 600 GeV")
    
    
    
    # line
    # line = ROOT.TLine(xmin+1,3.84,xmax-1,3.84)
#     line.SetLineColor(13)
#     line.SetLineWidth(1)
#     line.SetLineStyle(2)
#     line.Draw("same")
    
    
    #limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    #limits_file_300fb.write("%s & [%.2f,%.2f] \n"%(coupling,roots_2[0],roots_2[1]))
    
    if not os.path.isdir(args.InputDir + "/fitted_xsec_Marginalized"): os.mkdir(args.InputDir + "/fitted_xsec_Marginalized")

    c.SaveAs('%s/fitted_xsec_Marginalized/fit_xsec_Marginalized_%s.png'%(args.InputDir,coupling))
    c.SaveAs('%s/fitted_xsec_Marginalized/fit_xsec_Marginalized_%s.pdf'%(args.InputDir,coupling))
    
    print "done, saved in '%s/fitted_xsec__Marginalized/fit_xsec_Marginalized_%s.pdf"%(args.InputDir,coupling)
    
#limits_file.close()
#limits_file_300fb.close() 
    
    
    
    # f_, axarr_ = plt.subplots(2, sharex=True, figsize=(10,8))
#     
#     xmin = -20
#     xmax = +20
#     nsteps = 100
#     x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
#     #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
#     #y_centr = coeff_centr[0]+ coeff_centr[1]*x + coeff_centr[2]*x*x
#     y = func(x,*(coeff_centr))
#     y_up = func(x,*(coeff_centr + [i for i in errors_centr]))
#     y_down = func(x,*(coeff_centr - [i for i in errors_centr]))
#     axarr_[0].plot(x,y,c='r',label = 'fit')
#     axarr_[0].errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b",label='sample points')
#     axarr_[0].fill_between(x, y_down,y_up, facecolor='blue', alpha=0.3, label=r'fit $\pm 1\sigma$')
#     axarr_[0].fill_between(x, sm_xsec - sm_xsec_frac_error*sm_xsec, sm_xsec + sm_xsec_frac_error*sm_xsec, facecolor='green', alpha=0.3, label=r'SM 68\% CL')
#     axarr_[0].fill_between(x, sm_xsec - 2*sm_xsec_frac_error*sm_xsec, sm_xsec + 2*sm_xsec_frac_error*sm_xsec, facecolor='yellow', alpha=0.3, label=r'SM 95\% CL')
#     ymin = sm_xsec - 3*sm_xsec_frac_error*sm_xsec
#     ymax = sm_xsec + 6*sm_xsec_frac_error*sm_xsec
#     axarr_[0].axis([xmin, xmax, ymin, ymax])
#     axarr_[0].grid(True)
#     axarr_[0].set_xlabel(coupling, fontsize = 15)
#     axarr_[0].set_ylabel(r'cross section [pb]', fontsize = 15)
#     axarr_[0].legend(loc="upper right")
#     axarr_[0].text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff_centr[0],coeff_centr[1],coeff_centr[2]), fontsize=20, color="r")
#     
#     
#     coeff_centr, covmat_centr = curve_fit(func2,coupling_strengths,xsec,sigma=xsec_error)#absolute_sigma
#     errors_centr = sqrt(diag(covmat_centr))
#     
#     y_centr = func2(x,*(coeff_centr))
#     y_up = func2(x,*(coeff_centr + [i for i in errors_centr]))
#     y_down = func2(x,*(coeff_centr - [i for i in errors_centr]))
#     y_centr = [(i-sm_xsec)**2/((fracerr_measured_sm*sm_xsec)**2) for i in y_down]
#     ymin = -2
#     ymax = 10
#     axarr_[1].plot(x,y_centr,c="r",label=r'fitted $\chi^{2}$')
#     axarr_[1].axhline(y=3.84, linestyle="dashed", linewidth=2, color="navy", label=r'95\% CL')
#     p = P.fit(x, y_centr, 4)
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
#     
#     
#     
#     # y_centr = func2(x,*(coeff_centr))
# #     y_centr = [(i-sm_xsec)**2/((fracerr_measured_sm*sm_xsec)**2) for i in y_centr]
# #     ymax = 10#max(y_centr)
# #     ymin = 0
# #     #y_up = func4(x,*(coeff_up))
# #     #y_down = func4(x,*(coeff_down))
# #     #y_up = (coeff_centr[0]+errors_centr[0])*(1 + (coeff_centr[1]+errors_centr[1])*x + (coeff_centr[2]+errors_centr[2])*x*x)
# #     #y_down = (coeff_centr[0]-errors_centr[0])*(1 + (coeff_centr[1]-errors_centr[1])*x + (coeff_centr[2]-errors_centr[2])*x*x)
# #     #y_up = func4(x,*(coeff_centr + [i for i in errors_centr]))
# #     #y_down = func4(x,*(coeff_centr - [i for i in errors_centr]))
# #     
# #     plt.plot(x,y_centr,c="r",label='fit',linestyle="dotted")
# #     #plt.show()
# #     #plt.fill_between(x, y_down,y_up, facecolor='green', alpha=0.3, label=r'fitted 95% CL')
# #     plt.axhline(y=3.84, linestyle="dashed", linewidth=2, color="navy")
# #     p = P.fit(x, y_centr, 4)
# #     roots = sorted([i.real for i in (p - 3.84).roots() if i.imag == 0])
# #     print roots
# #     while len(roots)>2: roots = roots[1:-1]
# #     plt.axvline(x=roots[0], linestyle="dashed", linewidth=2, color="navy")
# #     plt.axvline(x=roots[1], linestyle="dashed", linewidth=2, color="navy")
# #     plt.legend(loc="upper right")
# # 
# # 
# #     # plot points with error bars
# #     #plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
# #     #if roots[0] < xmin: plt.axis([roots[0]-1, roots[1]+1, ymin, ymax])
# #     #else: 
# #     plt.axis([xmin, xmax, ymin, ymax])
# #     plt.grid(True)
# #     plt.xlabel(coupling, fontsize = 15)
# #     plt.ylabel('#chi^{2}/ndof', fontsize = 15)
# # 
# #     # draw some things on the canvas
# #     #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
# #     #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.5f) \ C^2 )$'%(coeff[0],coeff[1],coeff[2]), fontsize=12, color="r")
# #     #plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(sum(nevents.values())), fontsize=17, color="b")
# #     plt.text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
# #     plt.text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
# #     
#     #plt.savefig('%s/fit_ConfidenceLevelInterval_%s.png'%(args.OutputDir,c))
#     #plt.savefig('%s/fit_ConfidenceLevelInterval_%s.pdf'%(args.OutputDir,c))
#     
#     #plt.cla()
#     limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots[0],roots[1]))
#     
#     #print '%s: sigma=%.4f ( 1 + %fi  %.4fC + %.4f C^2 )$'%(coupling,coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0])
#     #print '%s: sigma=%.4f ( 1 + %.4fC + %.4f C^2 )$'%(coupling,coeff[0],coeff[1],coeff[2])
#     
#     plt.subplots_adjust(wspace=0.02, hspace=0)
#     
#     f_.savefig('%s/fit_xsec_%s.png'%(args.InputDir,coupling))
#     f_.savefig('%s/fit_xsec_%s.pdf'%(args.InputDir,coupling))
#     
#     f_.clf()
#     print "done, saved in %s/fit_xsec_%s.png"%(args.InputDir,coupling)
# 
# 
# limits_file.close()

# fsm_ = open(args.SMxsecDir+ "/cross_section_SM.txt", 'r')
# lines = fsm_.readlines()
# sm_xsec =  float(lines[-1].split(" ")[2])
# sm_xsec_error = float(lines[-1].split(" ")[3])
# fsm_.close()
# 
# f, axarr = plt.subplots(3, sharex=True, figsize=(12,12))
# x_ = range(len(coupling_names))
# plt.xticks(x_, coupling_names)
# axarr[0].errorbar(x_, p0, yerr = p0e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[0].set_ylabel("p0 [pb]", fontsize = 15)
# axarr[0].grid(True)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].set_ylim(sm_xsec - 3*sm_xsec_error, sm_xsec + 3*sm_xsec_error)
# textymin, textymax = axarr[0].get_ylim()
# axarr[0].text(x_[0]-1,textymax + (textymax - textymin)/float(10), r'$\sigma=p0 ( 1 + p1 \ C + p2 \ C^2 )$', fontsize=20)
# axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
# axarr[0].fill_between([-1] + x_ + [len(x_)+1], sm_xsec - sm_xsec_error, sm_xsec + sm_xsec_error, facecolor='red', alpha=0.3, label='SM effective cross section')
# axarr[0].legend(loc='best')
# axarr[1].errorbar(x_, p1, yerr = p1e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[1].set_ylabel("p1", fontsize = 15)
# axarr[1].grid(True)
# axarr[2].errorbar(x_, p2, yerr = p2e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# axarr[2].set_ylabel("p2", fontsize = 15)
# #axarr[2].set_xlabel("non-zero EFT coupling", fontsize = 15)
# axarr[2].grid(True)
# # axarr[3].errorbar(x_, p3, yerr = p3e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# # axarr[3].set_ylabel("p3", fontsize = 15)
# # #axarr[3].set_xlabel("non-zero EFT coupling", fontsize = 15)
# # axarr[3].grid(True)
# # axarr[4].errorbar(x_, p4, yerr = p4e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
# # axarr[4].set_ylabel("p4", fontsize = 15)
# # axarr[4].set_xlabel("non-zero EFT coupling", fontsize = 15)
# # axarr[4].grid(True)
# f.savefig('%s/fit_xsec_coeff_results.png'%(args.InputDir))
# f.savefig('%s/fit_xsec_coeff_results.pdf'%(args.InputDir))
    

