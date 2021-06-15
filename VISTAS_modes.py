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
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings

def LP_modes(wLength, nc, dn, rCore, nrho, rho, alpha):

    # STEP 1. calculate normalized cutoff frequencies and derive corresponding modes LPlm
    
    V = 2 * np.pi * rCore / wLength * nc * np.sqrt(2 * dn)    # normalized frequency for the specified waveguide

    lmax = 100
    mmax = 50
    ind = np.indices((lmax, mmax))
    lvec = ind[0].flatten('C')      # successive indices l and m flattened into vectors
    mvec = ind[1].flatten('C') + 1  # +1 to shift m (0 -> 1)

    VC = np.zeros((lmax, mmax))       # matrix with Vcutoff values stored along dimensions l and m (NB. l starts at zero, m at one)
    VC[0, 1 : mmax + 1] = np.round(sp.special.jn_zeros(1, mmax - 1), 5) # for l=0, Vcutoff as roots of J1(V) = 0, 1:max(mvec)-1 to account for added zero root not outputed by jn_zeros
    for l in range(1, lmax):
        VC[l, :] = np.round(sp.special.jn_zeros(l - 1, mmax), 5)        # for l!=0, Vcutoff given by the roots of Jl-1(V) = 0
    
    vc = VC.flatten('C')    # cutoff matrix flattened into a vector
    
    A = np.row_stack((lvec.astype(int), mvec.astype(int), vc))  # matrix with l-indices as first row, m-indices as second row, and Vcutoff as third row
    A = A[:, np.lexsort((A[1, :], A[2, :]))]                    # sort first by Vcutoff, then (for same Vcutoff) by m
    A = np.delete(A, np.where(A[2,:] > V), axis=1)              # remove columns where Vcutoff > V

    nm = A.shape[1]     # number of modes

    #LPlm = np.empty([nm], dtype = object)   # vector containing modes names as string elements
    LPlm = []   # list
    for m in range(nm):
        #print(f"LP{A[0,m].astype(int)}{A[1,m].astype(int)}: {np.round(A[2,m],3)}")
        # LPlm[m] = (f"LP{A[0,m].astype(int)}{A[1,m].astype(int)}")
        LPlm.append(f"LP{A[0,m].astype(int)}{A[1,m].astype(int)}")

    lvec = A[0, :].astype(int)
    mvec = A[1, :].astype(int)
      

    # STEP 2. Find roots u of dispersion equation, and format them into a 2D array along the axes l and m

    eigenValLP = lambda u, l: u * sp.special.jv(l+1, u) / sp.special.jv(l, u) - np.sqrt(V**2 - u**2) * sp.special.kv(
            l+1, np.sqrt(V**2 - u**2)) / sp.special.kv(l, np.sqrt(V**2 - u**2))  # Eigenvalue equation for LP modes
    
    u = np.zeros((max(lvec) + 1, max(mvec) + 1))                # 2D (l, m) array storing the roots of the eigenvalue equation for each LPlm mode
    warnings.filterwarnings('ignore', category=RuntimeWarning)  # prevents warning about np.sqrt(V**2 - u**2) when V < u

    for l in range(max(lvec) + 1):  # solve the eigenvalue equation for each of the l indices
        
        ul = np.array([0])  # roots for each l-index (rows of u)

        # run fsolve with multiple starting estimates between 1 and V (NB. the solution u must be < V)
        for u0 in np.linspace(1, V.astype(int) + 2, (V.astype(int) + 2) * alpha):  # the parameter 10 multiplies the number of starting point, thus improving convergence
            root, infodict, ier, msg = fsolve(lambda u: eigenValLP(u, l), u0, full_output=True)
            if ier == 1:
                ul = np.append(ul, np.round(root, 5)) # rounding enables subsequent suppression of duplicates using unique function

        ul = np.delete(ul, 0)       # remove roots = 0
        ul = np.delete(ul, ul <0)   # remove roots < 0
        ul = np.delete(ul, ul > V)  # remove roots > V
        ul = np.sort(np.unique(ul)) # remove duplicates and sort roots
        u[l, 1:len(ul) + 1] = ul

        # #problem visualization (plot of dispersion equation), useful for debugging and optimization
        # upts = np.linspace(0, V.astype(int) + 2, 100)
        # plt.plot(upts, eigenValLP(upts, l))
        # plt.ylim([-10, 10])
        # plt.grid()
        # plt.show()
        
    warnings.filterwarnings('default', category=RuntimeWarning) # switches warnings back on


    # STEP 3. Construct radial modal profiles based on roots found in step 2

    ur = np.zeros((nm, nrho))                 # radial modal intensity
    rBoundary = np.argmax(rho >= rCore)

    for m in range(nm):
        um = u[lvec[m], mvec[m]]
        if um != 0:
            wm = np.sqrt(V**2 - um**2)
            ur[m, :rBoundary] = sp.special.jv(lvec[m], um * rho[:, :rBoundary] / rCore) / sp.special.jv(lvec[m], um)   # profile in core
            ur[m, rBoundary:] = sp.special.kv(lvec[m], wm * rho[:, rBoundary:] / rCore) / sp.special.kv(lvec[m], wm)   # profile in cladding
        else:
            print(f'no root found for mode LP{lvec[m]}{mvec[m]}')

    return nm, lvec, LPlm, ur


def main():
    
    wLength = 858e-9
    rCore = 8e-6
    nc = 3.6
    dn = .001
    nrho = 100
    rho = np.linspace(0, 2.5*rCore, nrho)
    rho = rho[:, np.newaxis]    # np rank one array -> (1, nrho) vector

    nm, lvec, LPlm, ur = LP_modes(wLength, nc, dn, rCore, nrho, rho.T, alpha = 10)

    plt.plot(rho, ur.T) #plot field profiles
    plt.legend(LPlm)
    plt.show()


if __name__ == "__main__":
    main()