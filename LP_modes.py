##################################################################################
# MIT License
#
# Copyright (c) 2021 Marc Jungo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##################################################################################

import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import fsolve
import math

def LP_modes(wLength, nc, delta_n, r_act, r_cav, nrho, rho):
    
    import warnings

    ### Step 0. Define key parameters

    # # hard-coded parameters for the 20 lowest order modes
    lvec = np.array([0, 1, 2, 0, 3, 1, 4, 2, 0, 5, 3, 1, 6, 4, 7, 2, 0, 5, 8, 3])
    mvec = np.array([1, 1, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 1, 3, 4, 2, 1, 3])
    LPlm = np.array(['LP01', 'LP11', 'LP21', 'LP02', 'LP31', 'LP12', 'LP41', 'LP22', 'LP03', 'LP51', 'LP32', 'LP13', 'LP61', 'LP42', 'LP71', 'LP23', 'LP04', 'LP52', 'LP81', 'LP33'])
    Vcutoff = np.array([2.405, 3.832, 3.832, 5.136, 5.520, 6.380, 7.016, 7.016, 7.588,
                       8.417, 8.654, 8.771, 9.761, 9.936, 10.173, 10.173, 11.065, 11.086, 11.62]) #cutoff normalized frequencies for each of the successive LP modes

    # # code to calculate the normalized cutoff frequencies (hard-coded above)
    # VCutoff=np.zeros((max(lvec) + 1, max(mvec)))  # matrix with Vcutoff values stored along dimensions l and m (NB. l starts at zero, m at one)
    # VCutoff[0, 1:max(mvec)+1] = np.round(sp.special.jn_zeros(1, max(mvec)-1), 3)    # for l=0, Vcutoff as roots of J1(V) = 0, 1:max(mvec)-1 to account for added zero root not outputed by jn_zeros
    # for l in range(1, max(lvec)+1):
    #     VCutoff[l, :] = np.round(sp.special.jn_zeros(l-1, max(mvec)), 3)  # for l!=0, Vcutoff given by the roots of Jl-1(V) = 0
    # for m in range(len(lvec)):
    #     print("LP", lvec[m], mvec[m],": ", VCutoff[lvec[m], mvec[m]-1])

    V = 2 * np.pi * r_act / wLength * nc * np.sqrt(2 * delta_n)   #normalized frequency

    if V > max(Vcutoff):
        V = max(Vcutoff)
        nm=len(lvec)
        print('modes above LP33 (20th mode) are not calculated')
    else:
        nm = np.sum(Vcutoff <= V) + 1   #there are two modes for V>Vcutoff[0]=2.405, 3 for V>Vcutoff[1], etc.
    print(f'V-parameter: {round(V,3)}; # modes: {nm}')

    ### Step 1. Find roots u of dispersion equation, and format them into a 2D array along the axes l and m

    u=np.zeros((max(lvec) + 1, max(mvec) + 1))   #starts at zero
    warnings.filterwarnings('ignore', category=RuntimeWarning) #prevents warning about np.sqrt(V**2 - u**2) when V < u
    for l in range(max(lvec) + 1):              #solve the eigenvalue equation for each of the l indices

        # def eigenValLP(u): return u * sp.jvp(l, u, 1) / sp.jv(l, u) - np.sqrt(V**2 - u**2) * sp.kvp(
        # l, np.sqrt(V**2 - u**2), 1) / sp.kv(l + 1, np.sqrt(V**2 - u**2))  # Eigenvalue equation for LP modes

        def eigenValLP(u): return u * sp.jv(l+1, u) / sp.jv(l, u) - np.sqrt(V**2 - u**2) * sp.kv(
            l+1, np.sqrt(V**2 - u**2)) / sp.    kv(l, np.sqrt(V**2 - u**2))  # Eigenvalue equation for LP modes...
        
        u_tmp = np.array([0])  # roots for each l-index (rows)

        for u0 in range(1, math.ceil(V), 1):    #runs fsolve with multiple starting estimates between 1 and V (NB. the solution u must be <V)
            root, info, ier, mesg = fsolve(eigenValLP, u0, full_output=True)
            if ier == 1:
                u_tmp = np.append(u_tmp, np.round(root, 6)) #rounding enables suppression of duplicates using unique function

        u_tmp = np.sort(np.unique(u_tmp))       #remove duplicates and sort roots
        u_tmp = np.delete(u_tmp, 0)             #remove zeros
        u_tmp = np.delete(u_tmp, u_tmp > V)  # remove roots > V
        u[l, 1:len(u_tmp) + 1]=u_tmp

        # #problem visualization (plot of dispersion equation), useful for debugging and optimization
        # upts = np.linspace(0, math.ceil(V), 100)
        # plt.plot(upts, eigenValLP(upts))
        # plt.ylim([-10, 10])
        # plt.grid()
        # plt.show()

    # Step 2. Construct radial modal profiles

    ur = np.zeros((nm, nrho))                   #radial modal intensity
    ones_act = (rho <= r_act) * 1               #1 in active layer, 0 in cladding (* 1 to convert boot to int)
    ones_clad = (rho <= r_cav) * 1 - ones_act   #0 in active layer, 1 in cladding (* 1 to convert boot to int)
    for m in range(nm):
        um=u[lvec[m], mvec[m]]                  #2D array containing all roots, indexed according to the LPlm l and m indices
        wm=np.sqrt(V**2 - um**2)
        ur[m,:] = sp.jv(lvec[m], um * rho / r_act) / sp.jv(lvec[m], um) * ones_act + sp.kv(lvec[m], wm * rho / r_act) / sp.kv(lvec[m], wm) * ones_clad
    ur[:, 0] = ur[:, 1] #prevents edge artifacts
    
    plt.figure(1)
    plt.plot(rho.T, ur.T) #plot field profiles
    plt.legend(LPlm)
    plt.show()
    
    warnings.filterwarnings('default', category = RuntimeWarning)

    return(nm, lvec[0:nm], LPlm[0:nm], ur)

wLength = 858e-9
r_core = 4.5e-6
r_tot = 2 * r_core
nc = 3.6
dn = .001
nrho = 100
rho = np.linspace(0, r_tot, nrho)
rho = rho[:, np.newaxis]    # np rank one array -> (1, nrho) vector
nm, lvec, LPlm, ur = LP_modes(wLength, nc, dn, r_core, r_tot, nrho, rho.T)