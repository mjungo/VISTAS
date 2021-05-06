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

import matplotlib.pyplot as plt
import numpy as np
import time

from matplotlib import cm
from numba import njit
from scipy.optimize import fsolve
from scipy.optimize import check_grad
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.special import jn
from scipy.special import jn_zeros

import VISTASmodels

from LP_modes import LP_modes
from VCSEL_params import params

def main():
    
    # 1. key parameters --------------------------------------------------------------------------------
 
    (wl0, Rt, Rb, alpha_int, beta, Gamma, \
            tau_N, tau_esc , tau_cap, eta_i, rs, DN, \
                Ntr, gln, epsilon, \
                    r_act, r_cav, nc, ng, delta_n, L, V_cav, db, \
                        nrho, nphi) \
                            = params()

    q = 1.602e-19       # elementary charge [C]
    h = 6.626e-34       # Planck constant [Js]
    kB = 1.381e-23      # Boltzmann constant [J/K]
    c0 = 2.998e8        # vacuum speed of light [m/s]

    hvl = h * c0 * 1e9  # photon energy [J*m]
    vg = 100 * c0 / ng  # group velocity [cm/s]

    rho = np.linspace(0, r_cav, nrho)
    rho = rho[:, np.newaxis]    # np rank one array -> (1, nrho) vector
    phi = np.linspace(0, 2*np.pi, nphi)
    phi = phi[:, np.newaxis]    # np rank one array -> (1, nphi) vector

 
    # 2. normalized modes intensity profiles -----------------------------------------------------------

    # mode amplitude profiles
    nm, lvec, LPlm, ur = LP_modes(wl0 * 1e-9, nc, delta_n, r_act * 1e-2, r_cav * 1e-2, nrho, rho.T * 1e-2, alpha = 10)

    # normalized intensity profiles
    Ur = np.square(ur)  # field amplitude -> intensity
    norm = np.trapz(Ur * rho.T, rho.T) * 2 / r_cav**2
    norm = norm[:, np.newaxis]
    Ur = Ur / norm      # normalized intensity profile

    print()
    if str(input('Plot 2D mode profiles (y/n)? ')) == 'y':
        plot2D(Ur, LPlm, lvec, nm, rho, nrho, phi, nphi, nfig = 1)

    # size correction of parameters defined individually for each mode
    Rt = Rt[0:nm, :]
    Rb = Rb[0:nm, :]
    alpha_int = alpha_int[0:nm, :]
    epsilon = epsilon[0:nm, :]
    beta = beta[0:nm, :]
    Gamma = Gamma[0:nm, :]

    
    # 3. radial injection current profile --------------------------------------------------------------

    iact = (rho <= r_act) * 1   # "active" area, typically delimited by the oxide aperture
    ispread = (rho > r_act) * np.exp(-(rho - r_act) / rs)
    irho = iact + ispread
    norm = np.trapz(irho * rho, rho, axis = 0) * 2 / r_cav**2
    irho = irho / norm
    # plt.plot(rho, irho)
    # plt.show() 
    
    
    # 4. current and time characteristics --------------------------------------------------------------

    modformat = int(input('Step (0), single pulse (1), or large signal (2): '))
    tmax = float(input('  simulated time (ns): ')) * 1e-9
    tspan = [0, tmax]
    dt = 10e-12     # "reporting" timesteps (at which the solve_ivp is interpolated)
    dtFD = 1e-12    # timesteps for the finite difference evaluation
    teval = np.arange(0, tmax, dt)
    tevalFD = np.arange(0, tmax, dtFD)
    ct=len(teval)
    ctFD=len(tevalFD)

    It = np.ones((ct))      # used with solve_ivp solver (t_eval with typically 10ps timesteps)
    itFD = np.ones((ctFD))  # used with Finite Differences (FD) solver (t_eval with typically 1ps timesteps)

    if modformat == 0:
        Ion = float(input('  DC current (mA): ')) * 1e-3
        It = It * Ion
        itFD = itFD * Ion

    elif modformat == 1:
        Ion = float(input('  "on" current (mA): ')) * 1e-3
        Ioff = float(input('  "off" current (mA): ')) * 1e-3
        It = It * Ion
        It[int(ct/2):] = Ioff
        itFD = itFD * Ion
        itFD[int(ctFD/2):] = Ioff

    elif modformat == 2:
        Ion = float(input('  "on" current (mA): ')) * 1e-3
        Ioff = float(input('  "off" current (mA): ')) * 1e-3
        tb=float(input('  single bit duration (ns): ')) * 1e-9
        nb = int(tmax / tb)    # number of bits
        sequence = np.random.randint(2, size = nb)
        cb = int(tb / dt)                # number of time steps for one bit
        cbFD = int(tb / dtFD)
        It = It * Ion
        itFD = itFD * Ion
        for i in range(nb):
            It[i * cb : (i + 1) * cb] = It[i * cb : (i + 1) * cb] - (Ion - Ioff) * sequence[i]
            itFD[i * cbFD : (i + 1) * cbFD] = itFD[i * cbFD : (i + 1) * cbFD] - (Ion - Ioff) * sequence[i]


    # 5. coefficient arrays ----------------------------------------------------------------------------

    ni = int(input('Radial resolution (typically 7-20 terms): '))

    # carrier radial series expansion terms J0i = J0(gammai * rho / r_cav)
    gammai = np.asarray(jn_zeros(1, ni-1))   # scipy function that calculates the roots of the Bessel function
    gammai = np.insert(gammai, 0, 0, axis = 0)    # prepend 0 root, not included by jn_zero
    gammai = gammai[:, np.newaxis]
    
    J0i = jn(0, gammai * rho.T / r_cav)    # bessel exansion terms J0i (ni of them)

    # coefficients array for calculating the average carrier density Na over the active area Na
    c_act = np.zeros((1, ni))
    c_act = 2 / r_act**2 * np.trapz(J0i[:, 0:round(nrho * r_act / r_cav)].T * rho[0:round(nrho * r_act / r_cav), :], rho[0:round(nrho * r_act / r_cav), :], axis = 0)    # integrates in the active area only
    c_act = c_act[:, np.newaxis].T  # (1, ni)

    # coefficients array that captures the spatial overlap between current and carrier profiles
    c_inj = np.zeros((ni, 1))
    c_inj = np.trapz(J0i.T * irho * rho, rho, axis = 0)                 # (ni, )
    c_inj = c_inj[:, np.newaxis]                                        # (ni, 1)
    c_inj = c_inj * 2 * eta_i / q / V_cav / (r_cav * jn(0, gammai))**2  # (ni, 1)

    # coefficients array that describes non-stimulated carrier recombination (= 1 / tau_N)
    c_nst = np.ones((1, 1))
    c_nst = c_nst / tau_N

    # coefficients array that captures the overlap between carrier expansion terms to model diffusion
    c_diff = np.zeros((ni, 1))
    c_diff = DN * (gammai / r_cav)**2     # (ni, 1)

    # coefficients array that captures the overlap between carrier expansion terms and mode profiles
    c_st = np.zeros((nm, ni, ni))
    for m in range(nm):
        for i in range(ni):
            for j in range(ni):
                    prodTmp = J0i[j, :] * Ur[m, :] * J0i[i, :] * rho.reshape((-1,)) # reshape -> dim (nrho, )
                    prodTmp = prodTmp[:, np.newaxis]
                    c_st[m, i, j] = 2 * vg / (r_cav * jn(0, gammai[i]))**2 * np.trapz(prodTmp, rho, axis = 0)    

    # coefficients array that describes optical losses and outcoupling (= 1 / tau_S)
    c_sp = np.ones((nm, 1)) # assymettry (e.g. different internal losses) may be introduced for each mode
    alpha_m = 1 / L * np.log(1 / np.sqrt(Rt * Rb));         # mirror losses
    c_sp = c_sp * vg * (alpha_m + alpha_int) 
    
    F = (1 - Rt) / ((1 - Rt) + np.sqrt(Rt / Rb) * (1 - Rb)) # fraction of power emitted at the top facet
    eta_opt = F * alpha_m / (alpha_int + alpha_m)           # optical efficiency (eta_opt*eta_i=eta_d)  
    S2P = V_cav * hvl / wl0 * eta_opt * c_sp / Gamma * 1e3 # convert photon density to optical power (nm, 1)


    # 6. steady-state (LI) solution --------------------------------------------------------------------

    Imax = 10e-3                                # current range for LI characteristic
    dI = 0.1e-3                                 # current steps (must be small enough to ensure convergence)
    Icw = np.arange(0, Imax, dI)                # current vector
    NScw = np.zeros((ni + nm, Icw.shape[0]))    # NS solution matrix (ni+nm, Imax / dI)
    NScw[0:ni, 1] = (c_inj * Icw[1] / (c_nst + c_diff)).reshape((-1,))  # N-analytical solution (Rst = 0) for 2nd current step
    NScw[ni:ni + nm, 1] = np.maximum((Gamma * beta * c_nst * NScw[0, 1] / c_sp).reshape((-1,)), 0) # S-analytical solution (Rst = 0) for 2nd current step

    tcwSolverStart = time.time()    
    for i in range(2, Icw.shape[0], 1): # sol@0: 0, sol@1: NScw[:, 1] -> calculated analytically above
        args = (ni, nm, c_act, c_inj, Icw[i], c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp)
        NScw[:, i] = fsolve(VISTASmodels.cw_1D, NScw[:, i - 1], args = args, fprime = VISTASmodels.Jac_cw_1D)
    tcwSolverEnd = time.time()

    print()
    print(f'LI calculation: {np.round(tcwSolverEnd - tcwSolverStart, 3)}s')
    plotPower(Icw * 1e3, NScw[ni:,:], S2P, LPlm, xlabel = 'current (mA)')

    
    # 7a.solution of system of ODEs using 'solve_ivp' --------------------------------------------------

    NS = np.zeros((ni + nm, 1))
    NSinit = np.zeros((ni + nm))

    It = interp1d(x = teval, y = It, fill_value = "extrapolate") # current changes over time and must match the solve_ivp integration points
    args = (ni, nm, c_act, c_inj, It, c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp)
    
    tSolverStart = time.time()
    sol = solve_ivp(VISTASmodels.solver_1D, (0, tmax), NSinit, t_eval = teval, method = 'RK45', dense_output = True, vectorized = True, args = args)
    #sol = odeint(VISTASmodels.solver_1D, NSinit, t = teval, args = args, tfirst = True, Dfun = VISTASmodels.Jac_cw_1D)
    tSolverEnd = time.time()

    print()
    print(f'Solve_ivp main loop: {np.round(tSolverEnd - tSolverStart, 3)}s')
    plotPower(sol.t * 1e9, sol.y[ni:ni+nm], S2P, LPlm, xlabel = 'time (ns)')


    # 7b. solution of system of ODEs using finite differences ------------------------------------------

    # NFD = np.zeros((ni, ctFD))  # N0, N1, ..., NnN
    # Na = np.zeros((1, ctFD))    # average carrier density across the active area
    # SFD = np.zeros((nm, ctFD))  # multimode photon number matrix

    # Nto = np.zeros((ni, 1))
    # Sto = np.zeros((nm, 1))
    
    # tFDStart = time.time()
    
    # for ti in range(ctFD - 1):
    #     Nto[:, 0] = NFD[:, ti]
    #     Sto[:, 0] = SFD[:, ti]

    #     dNdt, dSdt = VISTASmodels.FD_1D(Nto, Sto, ni, nm, c_act, c_inj, itFD[ti], c_diff, gln, c_st, Ntr, epsilon, Gamma, beta, c_nst, c_sp) # rhs of system of ODEs

    #     Ntn = Nto + dtFD * dNdt             # Finite Differences step
    #     NFD[:, ti + 1] = Ntn.reshape((-1,))
        
    #     Stn = Sto + dtFD * dSdt             # Finite Differences step
    #     SFD[:, ti + 1] = Stn.reshape((-1,))

    # tFDEnd = time.time()
    
    # # subsample (smaller dtFD needed to ensure convergence, but too dense for post-processing)
    # NFD = NFD[:, ::int(dt/dtFD)]        
    # SFD = SFD[:, ::int(dt/dtFD)]
    # tevalFD = tevalFD[::int(dt/dtFD)]   # tevalFD = teval

    # print()
    # print(f'FD solution main loop: {np.round(tFDEnd - tFDStart, 3)}s')
    # plotPower(tevalFD * 1e9, SFD, S2P, LPlm, xlabel = 'time (ns)')


def plot2D(Ur, LPlm, lvec, nm, rho, nrho, phi, nphi, nfig):

    # create polar coordinates mesh
    P, R = np.meshgrid(phi, rho)

    # format 2D profiles
    UR = np.tile(Ur, (nphi, 1, 1))      # radial modal intensities (nphi, nm, nrho)
    UR = np.transpose(UR, (1, 2, 0))    # radial modal intensities (nm, nrho, nphi)
    UPc, UPs = np.zeros((nm, nrho, nphi)), np.zeros((nm, nrho, nphi))   # azimuthal modal intensities ("c" for cos(l*phi) modes, "s" for sin(l*phi) modes)
    for m in range(nm):
        UPc[m, :, :], UPs[m, :, :] = np.cos(lvec[m] * P)**2, np.sin(lvec[m] * P)**2
    Uc, Us = UPc * UR, UPs * UR

    # transform mesh into cartesian coordinates
    X, Y = R*np.cos(P), R*np.sin(P)

    # plot normalized intensity profiles
    for m in range(nm):
        fig = plt.figure(m+1)
        ax1 = plt.subplot(projection='3d')
        ax1.plot_surface(X*1e6, Y*1e6, Uc[m, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        #ax2 = plt.subplot()
        #ax2.pcolormesh(X*1e6, Y*1e6, Uc[m, :, :], antialiased=False, cmap=cm.seismic)
        #plt.axis('scaled')
        #fig_title = LPlm[m] + 'c, ' + 'normalized intensity (a.u.)'
        ax1.set_title(LPlm[m] + 'c' + '\nnormalized intensity (a.u.)')
        ax1.set_xlabel('cavity x (um)')
        ax1.set_ylabel('cavity y (um)')
        ax1.zaxis.set_visible(False)
        ax1.set_zticklabels([])
        plt.show()


def plotPower(x, S, S2P, modes, xlabel):
    P = S * S2P
    plt.plot(x, P.T)
    plt.legend(modes)
    plt.plot(x, np.sum(P, 0), 'k--')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(xlabel)
    plt.ylabel('optical output power (mW)')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()