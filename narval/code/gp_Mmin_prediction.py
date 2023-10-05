import numpy as np
import scipy as sp
import george
from george.kernels import *
import sys
sys.path.append("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax")
sys.path.insert(0, '/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASS_emulators/tools')
from gp_training import Invdisttree
from gp_training import *
from time import time
import scipy.interpolate

x = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax/GP_Mmin/HOD_design_np11_n5000_new_f_env_Msat_sigmalogM_fmax.dat")

x[:,0] = np.log10(x[:,0])
x[:,2] = np.log10(x[:,2])

x = x[:,0:-2]
xc = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASS_emulators/cosmology_camb_full.dat")

NewKernel = False
HODLR = False

Nsize1 = 0
Nsize2 = 50

N_hod_up = 0 # the parameter to change HODs: e.g. 0~4000 -> 1000~5000 when this number is 1000  

HH = np.array(range(0,4000))
HH  = HH.reshape(40, 100)
HH = HH + N_hod_up
HH = HH[:,Nsize1:Nsize2]
CC = range(39)
rr = np.empty((HH.shape[1]*len(CC), x.shape[1]+xc.shape[1]))
YY = np.empty(HH.shape[1]*len(CC))

pp = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax/GP_Mmin/Mmin_9bins_pp_Nsize_"+str(Nsize1)+"_"+str(Nsize2)+"+N_hod_up"+str(N_hod_up)+".dat")

HHmask = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax/GP_Mmin/HOD_mask_monopole.dat")
##################   find the mean of the data  #############

s2 = 0
for CID in CC:
    d = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax/GP_Mmin/training_mocks/Mmin_cosmo_"+str(CID)+".dat")
    HH3 = np.array(range(0, 100))
    for HID in HH3[Nsize1:Nsize2]:
        HID = int(HID)

        YY[s2] = d[HID]
        s2 = s2+1

Ymean = np.mean(YY)

##################  found the mean of the data ################

#GP_error = np.loadtxt("/home/zhai/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASS_emulators/"+NGAL+"/GP_satellite_fraction/satellite_fraction_error.dat")

GP_err = 0.01

y2 = np.empty((len(rr)*1))
ss2 = 0
for j in range(1):
    DC = j
    Ym = Ymean
    ss = 0
    yerr = np.zeros((len(rr)))
    y = np.empty((len(rr)))
    for CID in CC:
        d = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASSLOWZ_Msat_z07_fmax/GP_Mmin/training_mocks/Mmin_cosmo_"+str(CID)+".dat")
        
        HH3 = np.array(range(0, 100))
        for HID in HH3[Nsize1:Nsize2]:
            HID = int(HID)

            hh = HH[CID][HID]
            rr[ss,0:7]=xc[CID, :]
            rr[ss,7:16]=x[hh, :]
    
            y[ss] = d[HID]/Ymean
            yerr[ss] = GP_err
            
            y2[ss2] = y[ss]
            ss = ss+1
            ss2 = ss2+1

######

    p0 = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    k1 = ExpSquaredKernel(p0, ndim=len(p0))
    k2 = Matern32Kernel(p0, ndim=len(p0))
    k3 = ConstantKernel(0.1, ndim=len(p0))
    k4 = WhiteKernel(0.1, ndim=len(p0))
    k5 = ConstantKernel(0.1, ndim=len(p0))

    if NewKernel == False:
        kernel = k1*k5+k2+k3+k4
    else:
        kernel = k2+k5

    ppt = pp

    if j == 0:
        if HODLR == True:
            gp0 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp0 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp0.compute(rr, yerr)

        gp0.kernel.vector = ppt
        gp0.compute(rr, yerr)

    if j == 1:
        if HODLR == True:
            gp1 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp1 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp1.compute(rr, yerr)
        
        gp1.kernel.vector = ppt
        gp1.compute(rr, yerr)

    if j == 2:
        if HODLR == True:
            gp2 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp2 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp2.compute(rr, yerr)
        
        gp2.kernel.vector = ppt
        gp2.compute(rr, yerr)

    if j == 3:
        if HODLR == True:
            gp3 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp3 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp3.compute(rr, yerr)

        gp3.kernel.vector = ppt
        gp3.compute(rr, yerr)

    if j == 4:
        if HODLR == True:
            gp4 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp4 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp4.compute(rr, yerr)
        
        gp4.kernel.vector = ppt
        gp4.compute(rr, yerr)

    if j == 5:
        if HODLR == True:
            gp5 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp5 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp5.compute(rr, yerr)
        
        gp5.kernel.vector = ppt
        gp5.compute(rr, yerr)

    if j == 6:
        if HODLR == True:
            gp6 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp6 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp6.compute(rr, yerr)
        
        gp6.kernel.vector = ppt
        gp6.compute(rr, yerr)

    if j == 7:
        if HODLR == True:
            gp7 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp7 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp7.compute(rr, yerr)
        
        gp7.kernel.vector = ppt
        gp7.compute(rr, yerr)

    if j == 8:
        if HODLR == True:
            gp8 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp8 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp8.compute(rr, yerr)
        
        gp8.kernel.vector = ppt
        gp8.compute(rr, yerr)


def Prediction(param):
    tc = np.atleast_2d(param)
    mu0, cov0 = gp0.predict(y2[len(CC)*HH.shape[1]*0:len(CC)*HH.shape[1]*1], tc)
    pre = mu0*Ymean
    return pre
