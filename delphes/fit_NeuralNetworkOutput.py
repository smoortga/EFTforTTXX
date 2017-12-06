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

parser = ArgumentParser()

parser.add_argument('--InputDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/MODELSCAN_ttcc_inclusive_perOrder_DiLepton_LHCTOPWG", help='path to the directory of all restricted processes')
parser.add_argument('--ValidationDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_LHCTOPWG", help='path to the directory of all restricted processes')
parser.add_argument('--SMxsecDir', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0_backup/MODEL_ttcc_inclusive_SMRestrictedMassAndCabibbo", help='path to the directory of SM-like xsec measurement')
parser.add_argument('--TrainingFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_DILEPTON/training_output/model_checkpoint_save.hdf5", help='path to the directory of SM-like xsec measurement')
parser.add_argument('--ScalerFile', default = "/user/smoortga/Analysis/MG5_aMC_v2_6_0/CONVERTED_DELPHES_DILEPTON/training_output/scaler.pkl", help='path to the directory of SM-like xsec measurement')

args = parser.parse_args()



def NN_validate(filename,class_number=1,cut=0.,original_n_events = 50000):
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
    discr = model.predict(X)[:,class_number]
    nEvents = len(discr)
    print float(len(discr)), float(len(filename)*original_n_events), 100*float(len(discr))/float(len(filename)*original_n_events),"%"
    discr = discr[discr >= cut]
    #print 100*float(len(discr))/float(nEvents)
    #print ""
    #print filename, float(len(discr)),"/",float(len(filename)*original_n_events),"%"
    return float(len(discr))/float(len(filename)*original_n_events)

def extract_coupling(coupling,val):
    """ extract the numeric value of a coupling in the directory name"""
    if "minus" in val: return -float(val.replace(coupling+"minus","").replace("p","."))
    else: return float(val.replace(coupling,"").replace("p","."))
    
def func(x, a, b, c, d, e):
        return a*(1 + b*x + c*x*x + d*x*x*x + e*x*x*x*x)    
    

process_dict = {}
subdirs = [d for d in os.listdir(args.InputDir) if os.path.isdir(args.InputDir+"/"+d) and not "fitted" in d]
for subdir in subdirs:
    process_dict[subdir]= [dd for dd in os.listdir(args.InputDir+"/"+subdir) if os.path.isdir(args.InputDir+"/"+subdir+"/"+dd) and not "fitted_xsec" in dd]


validation_files = {}
files = [f for f in os.listdir(args.ValidationDir) if ".root" in f]
couplings = set([f.split("_")[0] for f in os.listdir(args.ValidationDir) if ".root" in f])
orders = set([f.split("_")[1] for f in os.listdir(args.ValidationDir) if ".root" in f])
for c in couplings:
    validation_files[c] = {}
    for o in orders:
        validation_files[c][o] = [f for f in files if c in f and o in f]

# coefficients with errors: sigma = p0 ( 1 + p1*C + p2*c^2 )
p0 = []
p0e = []
p1 = []
p1e = []
p2 = []
p2e = []
p3 = []
p3e = []
p4 = []
p4e = []
coupling_names = []

# SM cross section
fsm_ = open(args.SMxsecDir+ "/cross_section_SM.txt", 'r')
lines = fsm_.readlines()
sm_xsec =  float(lines[-1].split(" ")[2])
sm_xsec_error = float(lines[-1].split(" ")[3])
fsm_.close()

nevents=[50000]

result_dict = {}

for coupling,orders in process_dict.iteritems():
    result_dict[coupling] = {}
    for o in sorted(orders):
        print "Processing %s %s"%(coupling,o.split("_")[1])
        ##################
        #
        # Get cross section
        #
        ##################
        if os.path.isfile(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt"):
            f_ = open(args.InputDir + "/" + coupling + "/" + o + "/cross_section.txt", 'r')
            lines = f_.readlines()
            xsec = float(lines[-1].split(" ")[2])
            xsec_error = float(lines[-1].split(" ")[3])
        else: continue
        #xsec = 0.01
        #xsec_error = 0.001
        
        ##################
        #
        # NN output
        #
        ##################
        
        filelist = validation_files[coupling+"1p0"][o.split("_")[1]]
        #filelist = validation_files[coupling][o.split("_")[1]]
        frac_passing_evts = NN_validate(filelist,class_number=1,cut=0.3,original_n_events = nevents[0])
        result_dict[coupling][o.split("_")[1]] = [xsec,xsec_error,frac_passing_evts]
        
if not os.path.isdir(args.ValidationDir + "/validation_output"): os.mkdir(args.ValidationDir + "/validation_output")
pickle.dump(result_dict,open(args.ValidationDir + "/validation_output/result_dict.pkl",'wb'))
#result_dict = pickle.load(open(args.ValidationDir + "/validation_output/result_dict.pkl",'r'))




# find out the mean of the SM contributions as an average value
sm_values = []
for coupling,orders in result_dict.iteritems():
    for o,results in orders.iteritems():
        if "Order0" in o: sm_values.append(results[0]*results[2])
sm_xsec = np.mean(np.asarray(sm_values))
sm_xsec_error = np.std(np.asarray(sm_values))
    
# sm_xsec = 0.23
# sm_xsec_error = 0.01

fracerr_measured_sm = 0.2

# define a file for the limits to be stored: 
limits_file = open('%s/limits_NN_class1_cut0p3.txt'%(args.InputDir),"w")
#limits_file = open('%s/limits_nocut.txt'%(args.InputDir),"w")
limits_file.write("fractional error on SM measurement: %i %% \n"%int(100*fracerr_measured_sm))

for coupling,orders in result_dict.iteritems():
    coeff = [0,0,0,0,0]
    errors = [0,0,0,0,0]
    for o,results in orders.iteritems():
        o_int = int(o[-1])
        coeff[o_int] = results[0]*results[2]
        errors[o_int] = results[1]*results[2]
    for o,results in orders.iteritems():
        o_int = int(o[-1])
        if o_int != 0:
            #errors[i] = xsec_error[i]/xsec_error[0]
            errors[o_int] = sqrt( (errors[o_int]/coeff[0])**2 + (errors[0]*coeff[o_int]/(coeff[0]*coeff[0]))**2 )
            coeff[o_int] = coeff[o_int]/coeff[0]

    coupling_names.append(coupling)
    #print coeff, errors
    p0.append(coeff[0])
    p0e.append(errors[0])
    # p1.append(coeff[1]/coeff[0])
#     p1e.append(sqrt( (errors[1]/coeff[0])**2 - (errors[0]*coeff[1]/(coeff[0]*coeff[0]))**2 ))
#     p2.append(coeff[2]/coeff[0])
#     p2e.append(sqrt( (errors[2]/coeff[0])**2 - (errors[0]*coeff[2]/(coeff[0]*coeff[0]))**2 ))
    p1.append(coeff[1])
    p1e.append(errors[1])
    p2.append(coeff[2])
    p2e.append(errors[2])
    p3.append(coeff[3])
    p3e.append(errors[3])
    p4.append(coeff[4])
    p4e.append(errors[4])
    
    #plot fitted function
    xmin = -10
    xmax = +10
    ymin = min(p0)-(max(p0)- min(p0))/5.
    ymax = max(p0)+(max(p0)- min(p0))/5.
    nsteps = 30
    x = np.arange(xmin, xmax+ float(xmax-xmin)/float(nsteps), float(xmax-xmin)/float(nsteps))
    #y = coeff[0] + coeff[1]*x + coeff[2]*x*x
    y = coeff[0]*(1 + coeff[1]*x + coeff[2]*x*x + coeff[3]*x*x*x + coeff[4]*x*x*x*x)
    yup = (coeff[0]+errors[0])*(1 + (coeff[1]+errors[1])*x + (coeff[2]+errors[2])*x*x + (coeff[3]+errors[3])*x*x*x + (coeff[4]+errors[4])*x*x*x*x)
    ydown = (coeff[0]-errors[0])*(1 + (coeff[1]-errors[1])*x + (coeff[2]-errors[2])*x*x + (coeff[3]-errors[3])*x*x*x + (coeff[4]-errors[4])*x*x*x*x)
    # plt.plot(x,yup,c="b")
#     plt.plot(x,ydown,c="b")
    ymin = min(ydown)-(max(yup)-min(ydown))/5.
    ymax = max(yup)+(max(yup)-min(ydown))/5.
    plt.fill_between(x,ydown,yup, facecolor='blue', alpha = 0.2, label=r'fit $\pm 1\sigma$')
    plt.plot(x,y,c="r",label='fit')
    plt.fill_between(x, sm_xsec - fracerr_measured_sm*sm_xsec, sm_xsec + fracerr_measured_sm*sm_xsec, facecolor='green', alpha=0.3, label=r'SM $\pm$ '+str(100*fracerr_measured_sm)+'%')
    plt.fill_between(x, sm_xsec - 2*fracerr_measured_sm*sm_xsec, sm_xsec + 2*fracerr_measured_sm*sm_xsec, facecolor='yellow', alpha=0.3, label=r'SM $\pm$ '+str(200*fracerr_measured_sm)+'%')
    # get limits
    p = P.fit(x, ydown, 4)
    roots = sorted([i.real for i in (p - (sm_xsec + 2*fracerr_measured_sm*sm_xsec)).roots() if i.imag == 0])
    while len(roots)>2: roots = roots[1:-1]
    plt.axvline(x=roots[0], linestyle="dashed", linewidth=2, color="navy")
    plt.axvline(x=roots[1], linestyle="dashed", linewidth=2, color="navy")
    plt.legend(loc="upper right")
    
    
    # plot points with error bars
    # plt.errorbar(coupling_strengths, xsec, yerr = xsec_error, fmt= 'o', color = "b")
    if roots[0] < xmin: plt.axis([roots[0]-1, roots[1]+1, ymin, ymax])
    else: plt.axis([xmin, xmax, ymin, ymax])
    plt.grid(True)
    plt.xlabel(coupling, fontsize = 15)
    plt.ylabel('effective cross section (after MVA cut) [pb]', fontsize = 15)
    
    # draw some things on the canvas
    #plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.4f) \ C^2 )$'%(coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0]), fontsize=20, color="r")
    plt.text(xmin, ymax+(ymax-ymin)/float(20), r'$\sigma=%.4f ( 1 + (%.4f) \ C + (%.5f) \ C^2 + (%.6f) \ C^3 + (%.6f) \ C^4 )$'%(coeff[0],coeff[1],coeff[2],coeff[3],coeff[4]), fontsize=12, color="r")
    plt.text(xmin + (xmax-xmin)/float(20), ymax-1.5*(ymax-ymin)/float(20), r'%i events'%(nevents[0]), fontsize=17, color="b")
    plt.text(roots[0] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[0], fontsize=12, color="navy", ha = 'right')
    plt.text(roots[1] , ymin+(ymax-ymin)/float(20), "%.2f "%roots[1], fontsize=12, color="navy", ha = 'right')
    
    limits_file.write("%s & [%.2f,%.2f] \n"%(coupling,roots[0],roots[1]))
    
    #print '%s: sigma=%.4f ( 1 + %.4fC + %.4f C^2 )$'%(coupling,coeff[0],coeff[1]/coeff[0],coeff[2]/coeff[0])
    print '%s: sigma=%.4f ( 1 + %.4fC + %.4f C^2 + %.4f C^3 + %.4f C^4 )$'%(coupling,coeff[0],coeff[1],coeff[2],coeff[3],coeff[4])
    
    
    plt.savefig('%s/fit_xsec_%s.png'%(args.InputDir,coupling))
    plt.savefig('%s/fit_xsec_%s.pdf'%(args.InputDir,coupling))
    
    plt.close()
    print "done, saved in %s/fit_xsec_%s.png"%(args.InputDir,coupling)


limits_file.close()

# fsm_ = open(args.SMxsecDir+ "/cross_section_SM.txt", 'r')
# lines = fsm_.readlines()
# sm_xsec =  float(lines[-1].split(" ")[2])
# sm_xsec_error = float(lines[-1].split(" ")[3])
# fsm_.close()
# 
f, axarr = plt.subplots(5, sharex=True, figsize=(12,12))
x_ = range(len(coupling_names))
plt.xticks(x_, coupling_names)
axarr[0].errorbar(x_, p0, yerr = p0e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
axarr[0].set_ylabel("p0 [pb]", fontsize = 15)
axarr[0].grid(True)
axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
axarr[0].set_ylim(sm_xsec - 3*sm_xsec_error, sm_xsec + 3*sm_xsec_error)
textymin, textymax = axarr[0].get_ylim()
axarr[0].text(x_[0]-1,textymax + (textymax - textymin)/float(10), r'$\sigma=p0 ( 1 + p1 \ C + p2 \ C^2 + p3 \ C^3 + p4 \ C^4 )$', fontsize=20)
axarr[0].set_xlim(x_[0]-1, x_[-1]+1)
axarr[0].fill_between([-1] + x_ + [len(x_)+1], sm_xsec - sm_xsec_error, sm_xsec + sm_xsec_error, facecolor='red', alpha=0.3, label='SM effective cross section')
axarr[0].legend(loc='best')
axarr[1].errorbar(x_, p1, yerr = p1e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
axarr[1].set_ylabel("p1", fontsize = 15)
axarr[1].grid(True)
axarr[2].errorbar(x_, p2, yerr = p2e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
axarr[2].set_ylabel("p2", fontsize = 15)
#axarr[2].set_xlabel("non-zero EFT coupling", fontsize = 15)
axarr[2].grid(True)
axarr[3].errorbar(x_, p3, yerr = p3e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
axarr[3].set_ylabel("p3", fontsize = 15)
#axarr[3].set_xlabel("non-zero EFT coupling", fontsize = 15)
axarr[3].grid(True)
axarr[4].errorbar(x_, p4, yerr = p4e,fmt= 'o', color = "b", markersize=10, elinewidth=2)
axarr[4].set_ylabel("p4", fontsize = 15)
axarr[4].set_xlabel("non-zero EFT coupling", fontsize = 15)
axarr[4].grid(True)
f.savefig('%s/fit_xsec_coeff_results.png'%(args.InputDir))
f.savefig('%s/fit_xsec_coeff_results.pdf'%(args.InputDir))
    

