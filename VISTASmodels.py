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

# system of equations to be solved for calculation of steady-state (LI) characteristic
def cw_1D(y, *args): 
    
    ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gam_z, beta, c_nst, c_olo = args

    N = y[0 : ni]
    N = N[:, np.newaxis]
    S = y[ni : ni + nm]
    S = S[:, np.newaxis]
    
    Na = max(np.matmul(c_act, N), 1)                        # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # fitted time domain logarithmic gain factor
    #print(f'Na = {Na}; g0 = {g0}')

    Inj = c_inj * It                                        # current injection term
    Rnst = c_nst * N                                        # non-stimulated carrier recombination
    Diff = c_diff * N                                       # carrier diffusion term
    Rst = g0 * (np.squeeze(np.matmul(c_st, N)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination term
    Rolo = c_olo * S                                        # optical losses and outcoupling term

    Nsol = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    Ssol = -Rolo + Gam_z * beta * Rnst[0, :] + Gam_z * Rst[:, 0][:, np.newaxis]
    NSsol = np.concatenate([Nsol, Ssol])
    NSsol = NSsol.reshape((-1,))

    return NSsol

# Jacobian (ni + nm, ni + nm) used by fsolve to calculate the steady-state (LI) characteristic
def Jac_cw_1D(y, *args):

    ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gam_z, beta, c_nst, c_olo = args

    N = y[0 : ni]
    N = N[:, np.newaxis]
    S = y[ni : ni + nm]
    S = S[:, np.newaxis]

    Na = max(np.matmul(c_act, N), 1)                        # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # fitted time domain logarithmic gain factor
    #print(f'Na = {Na}; g0 = {g0}')

    diag_ni = np.zeros((ni, ni), int)
    np.fill_diagonal(diag_ni, 1)
    diag_nm = np.zeros((nm, nm), int)
    np.fill_diagonal(diag_nm, 1)
    col = np.zeros((nm, ni))
    col[:, 0] = 1

    Ngain = g0 * (np.squeeze(np.matmul(c_st, N)) - c_st[:, :, 0] * Ntr)
    Scompr = 1 + epsilon * S
    Sterm = S / Scompr

    d11 = -(c_nst + c_diff) * diag_ni - g0 * np.sum(c_st * Sterm[:,:,np.newaxis] , 0)   # (ni,ni)
    d12 = -(Ngain / Scompr**2).T                                                        # (ni,nm)
    d21 = Gam_z * beta * c_nst * col + Gam_z * g0 * c_st[:, 0, :] * Sterm               # (nm,ni)
    d22 = (-c_olo + Gam_z * Ngain[:, 0] / Scompr**2) * diag_nm                          # (nm,nm)

    col1 = np.concatenate((d11,d21), 0)
    col2 = np.concatenate((d12,d22), 0)
    
    return np.concatenate((col1, col2), 1)


# system of ODEs to be solved for dynamic response calculation using solve_ivp solver
def solver_1D(t, y, *args):

    ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gam_z, beta, c_nst, c_olo = args 

    N = y[0 : ni, 0]
    N = N[:, np.newaxis]
    S = y[ni : ni + nm, 0]
    S = S[:, np.newaxis]

    Na = max(np.matmul(c_act, N), 1)                        # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # fitted time domain logarithmic gain factor

    Inj = c_inj * It(t)                                     # current injection
    Rnst = c_nst * N                                        # non-stimulated carrier recombination
    Diff = c_diff * N                                       # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, N)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination
    Rolo = c_olo * S                                        # optical losses and outcoupling term

    dNdt = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rolo + Gam_z * beta * Rnst[0, :] + Gam_z * Rst[:, 0][:, np.newaxis]
    
    return np.concatenate([dNdt, dSdt])


# system of ODEs to be solved for dynamic response calculation using finite differences
def FD_1D(Nto, Sto, ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gam_z, beta, c_nst, c_olo): 
    
    Na = max(np.matmul(c_act, Nto), 1)                        # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # fitted time domain logarithmic gain factor

    Inj = c_inj * It                                        # current injection
    Rnst = c_nst * Nto                                      # non-stimulated carrier recombination
    Diff = c_diff * Nto                                     # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, Nto)) \
         - c_st[:, :, 0] * Ntr) * Sto / (1 + epsilon * Sto) # stimulated recombination
    Rolo = c_olo * Sto                                      # optical losses and outcoupling term

    dNdt = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rolo + Gam_z * beta * Rnst[0, :] + Gam_z * Rst[:, 0][:, np.newaxis]
    
    return dNdt, dSdt