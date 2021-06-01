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

import json
import matplotlib.pyplot as plt
import numpy as np
import time

from matplotlib import cm
#from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
#from scipy.optimize import check_grad
#from scipy.signal import savgol_filter
from scipy.special import jn
from scipy.special import jn_zeros

import VISTASmodels

from LP_modes import LP_modes


def VISTAS1D(sp, vp):
    print()
    # 1. key parameters --------------------------------------------------------------------------------
    q = 1.602e-19               # elementary charge [C]
    h = 6.626e-34               # Planck constant [Js]
    #kB = 1.381e-23             # Boltzmann constant [J/K]
    c0 = 2.998e8                # vacuum speed of light [m/s]

    hvl = h * c0 * 1e9          # photon energy [J*m]
    vg = 100 * c0 / vp['ng']    # group velocity [cm/s]

    r_cav = 3 * vp['r_ox']                              # optical field and carriers extend beyond the active layer, delimited by the oxide
    V_cav = np.pi * r_cav**2 * vp['nqw'] * vp['dqw']    # cavity volume (cm-3)
    Gam_z = vp['dqw'] * vp['nqw'] / vp['Leff']          # longitudinal optical confinement factor 

    nrho = 100  # radial steps
    nphi = 200  # azimuthal steps

    rho = np.linspace(0, r_cav, nrho)
    rho = rho[:, np.newaxis]                # np rank one array -> (1, nrho) vector (rank 2 array)
    phi = np.linspace(0, 2*np.pi, nphi)
    phi = phi[:, np.newaxis]                # np rank one array -> (1, nphi) vector (rank 2 array)

 
    # 2. normalized modes intensity profiles -----------------------------------------------------------
    # mode amplitude profiles
    
    tStart = time.time()    
    nm, lvec, LPlm, ur = LP_modes(vp['wl0']* 1e-9, vp['nc'], vp['delta_n'], vp['r_ox'] * 1e-2, nrho, rho.T * 1e-2, alpha = 10)
    tEnd = time.time()
    print(f'Mode profiles computation: {np.round(tEnd - tStart, 3)}s')

    # normalized intensity profiles
    Ur = np.square(ur)  # field amplitude -> intensity
    norm = np.trapz(Ur * rho.T, rho.T) * 2 / r_cav**2
    norm = norm[:, np.newaxis]
    Ur = Ur / norm      # normalized intensity profile

    # size correction of parameters defined individually for each mode
    alpha_m = 1 / vp['Leff'] * np.log(1 / np.sqrt(vp['Rt'] * vp['Rb'])) * np.ones((nm, 1))   # mirror losses
    alpha_i = vp['alpha_i'] * np.ones((nm, 1))
    epsilon = vp['epsilon'] * np.ones((nm, 1))
    beta = vp['beta'] * np.ones((nm, 1))
    Gam_z = Gam_z * np.ones((nm, 1))

    
    # 3. radial injection current profile --------------------------------------------------------------
    iact = (rho <= vp['r_ox']) * 1   # "active" area, typically delimited by the oxide aperture
    ispread = (rho > vp['r_ox']) * np.exp(-(rho - vp['r_ox']) / vp['rs'])
    irho = iact + ispread
    norm = np.trapz(irho * rho, rho, axis = 0) * 2 / r_cav**2
    irho = irho / norm
    # plt.plot(rho, irho)
    # plt.show() 
    

    # 4. coefficient arrays ----------------------------------------------------------------------------
    tStart = time.time()     
    # carrier radial series expansion terms J0i = J0(gammai * rho / r_cav)
    gammai = np.asarray(jn_zeros(1, sp['ni']-1))      # scipy function that calculates the roots of the Bessel function
    gammai = np.insert(gammai, 0, 0, axis = 0)  # prepend 0 root, not included by jn_zero
    gammai = gammai[:, np.newaxis]
    J0i = jn(0, gammai * rho.T / r_cav)         # bessel expansion terms J0i (sp['ni'] of them)

    # coefficients array for calculating the average carrier density Na over the active area Na
    c_act = np.zeros((1, sp['ni']))
    c_act = 2 / vp['r_ox']**2 * np.trapz(J0i[:, 0:round(nrho * vp['r_ox'] / r_cav)].T * rho[0:round(nrho * vp['r_ox'] / r_cav), :], rho[0:round(nrho * vp['r_ox'] / r_cav), :], axis = 0)    # integrates in the active area only
    c_act = c_act[:, np.newaxis].T  # (1, sp['ni'])

    # coefficients array that captures the spatial overlap between current and carrier profiles
    c_inj = np.zeros((sp['ni'], 1))
    c_inj = np.trapz(J0i.T * irho * rho, rho, axis = 0)                 # (sp['ni'], )
    c_inj = c_inj[:, np.newaxis]                                        # (sp['ni'], 1)
    c_inj = c_inj * 2 * vp['eta_i'] / q / V_cav / (r_cav * jn(0, gammai))**2  # (sp['ni'], 1)

    # coefficients array that describes non-stimulated carrier recombination (= 1 / tau_N)
    c_nst = np.ones((1, 1))
    c_nst = c_nst / vp['tau_N']

    # coefficients array that captures the overlap between carrier expansion terms to model diffusion
    c_diff = np.zeros((sp['ni'], 1))
    c_diff = vp['DN'] * (gammai / r_cav)**2     # (sp['ni'], 1)

    # coefficients array that captures the overlap between carrier expansion terms and mode profiles
    c_st = np.zeros((nm, sp['ni'], sp['ni']))
    for m in range(nm):
        for i in range(sp['ni']):
            for j in range(sp['ni']):
                    prodTmp = J0i[j, :] * Ur[m, :] * J0i[i, :] * rho.reshape((-1,)) # reshape -> dim (nrho, )
                    prodTmp = prodTmp[:, np.newaxis]
                    c_st[m, i, j] = 2 * vp['Gam_r'] * vg / (r_cav * jn(0, gammai[i]))**2 * np.trapz(prodTmp, rho, axis = 0)    

    # coefficients array that describes optical losses and outcoupling (= 1 / tau_S)
    c_olo = np.ones((nm, 1)) # assymettry (e.g. different internal losses) may be introduced for each mode
    c_olo = c_olo * vg * (alpha_m + alpha_i) 
    
    # parameters for optical power calculation
    F = (1 - vp['Rt']) / ((1 - vp['Rt']) + np.sqrt(vp['Rt'] / vp['Rb']) * (1 - vp['Rb'])) # fraction of power emitted at the top facet
    eta_opt = F * alpha_m / (alpha_i + alpha_m)             # optical efficiency (eta_opt*eta_i=eta_d)  
    S2P = eta_opt * hvl / vp['wl0'] * V_cav * c_olo / Gam_z * 1e3 # conversion factor photon density -> optical power (nm, 1)
    tEnd = time.time()
    print(f'Overlap parameters computation: {np.round(tEnd - tStart, 3)}s')
    
    # 5. current and time characteristics --------------------------------------------------------------
    if sp['odeSolver'] == 'Finite Diff.':
        dt = sp['dtFD']
    else:
        dt = sp['dt']

    teval = np.arange(0, sp['tmax'], dt)
    ct = len(teval)

    if sp['modFormat'] == 'step':
        It = np.ones((ct)) * sp['Ion']      # current vector
        It[0] = 0   # to calculate NSinit@I(t0)=0

    elif sp['modFormat'] == 'pulse':
        It = np.ones((ct)) * sp['Ion']      # current vector
        It[int(ct/2):] =  sp['Ioff']
        It[0] = sp['Ioff']   # to compute NSinit@I(t0)=0

    elif sp['modFormat'] == 'random bits':
        nb = int(sp['tmax'] / sp['tb'])     # number of bits
        sequence = np.random.randint(2, size = nb)
        cb = int(sp['tb'] / dt)             # number of time steps for one bit
        It = np.ones((ct)) * sp['Ion']      # current vector
        for i in range(nb):
            It[i * cb : (i + 1) * cb] = It[i * cb : (i + 1) * cb] - (sp['Ion'] - sp['Ioff']) * sequence[i]

    elif sp['modFormat'] == 'small signal':
        ct = 2**14 + 1     # 8192 points -> df = 1/dt/ct = 1/tmax = 12 MHz (for dt = 1ps)
        sp['tmax'] = ct * dt
        teval = np.arange(0, sp['tmax'], dt)
        It = np.ones((ct)) * (sp['Ion'] + sp['Iss'])
        It[0] = sp['Ion']             # to compute NSinit@I(t0)=Ion

    # 6. steady-state (LI) solution --------------------------------------------------------------------
    tStart = time.time() 
    Imax = 10e-3                                # current range for LI characteristic
    dI = 0.1e-3                                 # current steps (must be small enough to ensure convergence)

    Icw = np.arange(0, Imax, dI)                # current vector
    NScw = np.zeros((sp['ni'] + nm, Icw.shape[0]))    # NS solution matrix (sp['ni']+nm, Imax / dI)
    NScw[0:sp['ni'], 1] = (c_inj * Icw[1] / (c_nst + c_diff)).reshape((-1,))  # N-analytical solution (Rst = 0) for 2nd current step
    NScw[sp['ni']:sp['ni'] + nm, 1] = np.maximum((Gam_z * beta * c_nst * NScw[0, 1] / c_olo).reshape((-1,)), 0) # S-analytical solution (Rst = 0) for 2nd current step
   
    for i in range(2, Icw.shape[0], 1): # sol@0: 0, sol@1: NScw[:, 1] -> calculated analytically above
        args = (sp['ni'], nm, c_act, c_inj, Icw[i], c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, Gam_z, beta, c_nst, c_olo)
        NScw[:, i] = fsolve(VISTASmodels.cw_1D, NScw[:, i - 1], args = args, fprime = VISTASmodels.Jac_cw_1D)

    NScwInterp = interp1d(x = Icw, y = NScw, fill_value = "extrapolate")    # extrapolates the LI characteristic -> NScwInterp(I)

    tEnd = time.time()
    print(f'LI solution: {np.round(tEnd - tStart, 3)}s')

    # 7. solver main loop ------------------------------------------------------------------------------
    tStart = time.time() 
    N = np.zeros((sp['ni'], ct))            # N0, N1, ..., NnN
    S = np.zeros((nm, ct))                  # multimode photon number matrix
    NSinit = NScwInterp(It[0])              # result of steady-state (LI) characteristic (extrapolated) as starting point to avoid initial fluctuations

    if sp['odeSolver'] == 'Finite Diff.':       # 7a. solution of system of ODEs using finite differences with constant time-step

        Nto = np.zeros((sp['ni'], 1))           # temp variable for passed to FD calculation function
        Sto = np.zeros((nm, 1))                 # temp variable for passed to FD calculation function

        N[:, 0] = NSinit[:sp['ni']]
        S[:, 0] = NSinit[sp['ni']:]
     
        for ti in range(ct - 1):
            Nto[:, 0] = N[:, ti]
            Sto[:, 0] = S[:, ti]

            dNdt, dSdt = VISTASmodels.FD_1D(Nto, Sto, sp['ni'], nm, c_act, c_inj, It[ti], c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, Gam_z, beta, c_nst, c_olo) # rhs of system of ODEs

            Ntn = Nto + dt * dNdt               # Finite Differences step
            N[:, ti + 1] = Ntn.reshape((-1,))
            
            Stn = Sto + dt * dSdt               # Finite Differences step
            S[:, ti + 1] = Stn.reshape((-1,))
        tEnd = time.time()
        print(f'FD ODE solution: {np.round(tEnd - tStart, 3)}s')

        tStart = time.time()
        if sp['modFormat'] == 'small signal':
            f, H = freqResp(S[:,1:], S2P, ct-1, dt)   # compute frequency response
            tEnd = time.time()
            print(f'Frequency response calculation: {np.round(tEnd - tStart, 3)}s')

        # subsampling: smaller timestep needed to ensure convergence, but too dense for post-processing  
        N = N[:, ::int(sp['dt']/dt)]
        S = S[:, ::int(sp['dt']/dt)]
        teval = teval[::int(sp['dt']/dt)]
        ct = len(teval)
        dt = sp['dt']
 
    else:   # 7b.solution of system of ODEs using 'solve_ivp'

        tStart = time.time()        
        It = interp1d(x = teval, y = It, fill_value = "extrapolate") # current changes over time and must match the solve_ivp integration points
        args = (sp['ni'], nm, c_act, c_inj, It, c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, Gam_z, beta, c_nst, c_olo)
        
        sol = solve_ivp(VISTASmodels.solver_1D, (0, sp['tmax']), NSinit, t_eval=teval, method=sp['odeSolver'], dense_output=True, vectorized=True, args=args, rtol=1e-6, atol=1e-6)
        #sol = odeint(VISTASmodels.solver_1D, NSinit, t = teval, args = args, tfirst = True, Dfun = VISTASmodels.Jac_cw_1D)

        N = sol.y[:sp['ni']]
        S = sol.y[sp['ni']:sp['ni']+nm]
        tEnd = time.time()
        print(f'solve_ivp ODE solution: {np.round(tEnd - tStart, 3)}s')

        tStart = time.time()
        if sp['modFormat'] == 'small signal':
            f, H = freqResp(S[:,1:], S2P, ct-1, dt)   # compute frequency response
            tEnd = time.time()
            print(f'Frequency response calculation: {np.round(tEnd - tStart, 3)}s')

    # 8. visualization ---------------------------------------------------------------------------------
    if sp['Uxyplot'] == 1:  # 2D mode profiles (cosine azimuthal distribution) Ur(x,y)
        plot2D(Ur, LPlm, lvec, nm, rho*1e-2, nrho, phi, nphi, nfig = 1)              
    
    if sp['PIplot'] == 1:   # steady-state LI characteristic Popt(I)
        plotPower(Icw * 1e3, NScw[sp['ni']:,:], S2P, LPlm, xlabel = 'current (mA)') 

    if sp['Ptplot'] == 1:   # dynamic response Popt(t)
        plotPower(teval * 1e9, S, S2P, LPlm, xlabel = 'time (ns)')
        
    if sp['modFormat'] == 'small signal' and sp['Hfplot'] == 1: # small signal response H(f)
        plotH(f, H)                                                                

    # 9. simulation results in dictionary sr -----------------------------------------------------
    sr={ 
            # "phi": phi.tolist(),    # dictionary can't store numpy arrays -> conversion to list
            # "nphi": nphi,
            # "rho": rho.tolist(),
            # "nrho": nrho,
            # "Icw": Icw.tolist(),
            # "NScw": NScw.tolist(),
            # "teval": teval.tolist(),
            # "S2P": S2P.tolist(),
            # "lvec": lvec.tolist(),
            # "f": f.tolist(),
            # # "H": H.tolist(),
            "nm": nm,
            "LPlm": LPlm,
            "ur": ur.tolist(),
            "N": N.tolist(),
            "S": S.tolist(),
        }
    print()
    return sr


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
    X, Y = R * np.cos(P), R * np.sin(P)

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
        #plt.show(block = False)
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


def freqResp(S, S2P, ct, dt):

    # 0. frequency vector (GHz)
    f = np.linspace(start=0, stop=int(1/2/dt), num=int(ct/2))*1e-9
    # 1. response to small signal step
    P = np.sum(S * S2P, 0)
    # 2. derivative of unit step = impulse: H(f) = Y(f) = FFT(yimp(t))) = FFT(d/dt ystep(T)
    dP = np.gradient(P, dt)    # dP = np.diff(P) / dt as alternative
    # 3. FFT to transform to frequency domain (-> complex)
    DP = np.fft.fft(dP)
    # 4. two-sided -> one-sided -> multiply amplitude of first half by 2, discard the rest
    DP = DP[:int(ct/2)] * 2
    #DP[0] = DP[0]/2
    # 5. compute amplitude as product of DP and its complex conjugate, normalize, remove 0 complex element
    H = DP * np.conj(DP) / ct
    H = H / H[0]
    H = H.real
    # 6. convert do dB
    H = 10 * np.log10(H)
    # (7. optional: filter to remove artefacts, e.g. in case rtol is too large in solve_ivp)
    #H = savgol_filter(H, window_length=(ct/2), polyorder=4, mode='mirror')

    return f, H


def plotH(f, H):

    plt.plot(f, H)
    plt.xlim(xmin=0, xmax=10)
    plt.ylim(ymin=-10, ymax=15)
    plt.xlabel("Frequency (GHz)")
    plt.ylabel('Frequency response (dB)')
    plt.grid()
    plt.show()


def main():
    
    with open('last_params.json', 'r') as openfile:
        params = json.load(openfile)    # dictionary
    sp, vp = params['simParams'], params['vcselParams']
    del params

    sr = VISTAS1D(sp, vp)


if __name__ == "__main__":
    main()