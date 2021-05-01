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

# system of equations to be solved for calculation of LI characteristic
def cw_1D(NS, ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp): 
    
    N = NS[0 : ni]
    N = N[:, np.newaxis]
    S = NS[ni : ni + nm]
    S = S[:, np.newaxis]
    
    #g0 = 5e-16
    Na = np.matmul(c_act, N)                                # average carrier density over the active area
    g0 = gln * np.log((Na + 1) / Ntr) / (Na-Ntr)            # fitted time domain logarithmic gain factor

    Inj = c_inj * It                                        # current injection term
    Diff = c_diff * N                                       # carrier diffusion term
    Rst = g0 * (np.squeeze(np.matmul(c_st, N)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination term
    Rnst = c_nst * N                                        # non-stimulated photon recombination
    Rsp = c_sp * S                                          # spontaneous photon recombination

    Nsol = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    Ssol = -Rsp + Gamma * beta * Rnst[0, :] + Gamma * Rst[:, 0][:, np.newaxis]
    NSsol = np.concatenate([Nsol, Ssol])
    NSsol = NSsol.reshape((-1,))

    return NSsol


# system of ODEs to be solved for dynamic response calculation using solve_ivp solver
def solver_1D(ti, NS, ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp): 
    
    N = NS[0 : ni, 0]
    N = N[:, np.newaxis]
    S = NS[ni : ni + nm, 0]
    S = S[:, np.newaxis]
    
    #g0 = 5e-16
    Na = np.matmul(c_act, N)                                # average carrier density over the active area
    g0 = gln * np.log((Na + 1) / Ntr) / (Na-Ntr)            # fitted time domain logarithmic gain factor

    Inj = c_inj * It(ti)                                    # current injection
    Diff = c_diff * N                                       # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, N)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination
    Rnst = c_nst * N                                        # non-stimulated photon recombination
    Rsp = c_sp * S                                          # spontaneous photon recombination

    dNdt = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rsp + Gamma * beta * Rnst[0, :] + Gamma * Rst[:, 0][:, np.newaxis]
    dNSdt = np.concatenate([dNdt, dSdt])

    return dNSdt


# system of ODEs to be solved for dynamic response calculation using finite differences
def FD_1D(Nto, Sto, ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp): 
    
    #g0 = 5e-16
    Na = np.matmul(c_act, Nto)                              # average carrier density over the active area
    g0 = gln * np.log((Na + 1) / Ntr) / (Na-Ntr)            # fitted time domain logarithmic gain factor

    Inj = c_inj * It                                        # current injection
    Diff = c_diff * Nto                                     # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, Nto)) \
         - c_st[:, :, 0] * Ntr) * Sto / (1 + epsilon * Sto) # stimulated recombination
    Rnst = c_nst * Nto                                      # non-stimulated photon recombination
    Rsp = c_sp * Sto                                        # spontaneous photon recombination
    Rsp = c_sp * S                                          # spontaneous photon recombination

    dNdt = Inj - Rnst - Diff - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rsp + Gamma * beta * Rnst[0, :] + Gamma * Rst[:, 0][:, np.newaxis]
    
    return dNdt, dSdt