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
# import matplotlib.pyplot as plt
import numpy as np
import time

#from scipy.integrate import odeint
from scipy.fft import fft, ifft, fftshift, fftfreq
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
# from scipy.optimize import check_grad
# from scipy.signal import savgol_filter
from scipy.special import jn
from scipy.special import jn_zeros

from VISTAS_modes import LP_modes


# system of equations to be solved for calculation of steady-state (LI) characteristic
def cw_1D_transp(y, *args): 
    
    SCHtransp, nNb, nNw, nS, Vr, c_act, c_injb, c_injw, It, c_diff, gln, c_st, Ntr, epsilon, GamZ, beta, tauNb, tauNw, tauCap, tauEsc, taub, tauw, tauS = args

    Nb=y[0]
    Nw = y[1:nNb+nNw]
    Nw = Nw[:, np.newaxis]
    S = y[nNb+nNw:]
    S = S[:, np.newaxis]
    
    Na = max(np.matmul(c_act, Nw), 1)                       # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # instantaneous linear gain coefficient, fitted to logarithmic gain function (-> time domain)
    #print(f'Na = {Na}; g0 = {g0}')

    Injb = c_injb * It                                      # current injection
    
    if SCHtransp == True:
        Cap = Nb / tauCap                                   # carrier capture term from the barriers into the QWs
        Esc = Nw / tauEsc                                   # carrier escape term from the QWs into the barriers
        Rnstb = Nb / tauNb                                  # non-stimulated carrier recombination in the barriers
    else:  
        Cap = Injb                                          # carrier capture term from the barriers into the QWs
        Esc = np.zeros((nNw, 1))                            # carrier escape term from the QWs into the barriers
        Rnstb = 0                                           # non-stimulated carrier recombination in the barriers

    Injw = c_injw / Vr * Cap                                # current injection into the wells
    Rnstw = Nw / tauNw                                      # non-stimulated carrier recombination
    Diff = c_diff * Nw                                      # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, Nw)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination
    Rsp = GamZ * beta * Rnstw[0, :]                         # spontaneous emission (spread over modes based on their relative power)
    Rolo = S / tauS                                         # optical losses and outcoupling term

    Nbsol = Injb - Rnstb - Cap + Vr * Esc[0,:][:, np.newaxis]
    Nwsol = Injw - Rnstw - Diff - Esc - np.sum(Rst, 0)[:, np.newaxis]
    Ssol = -Rolo + Rsp + GamZ * Rst[:, 0][:, np.newaxis]
    NSsol = np.concatenate([Nbsol, Nwsol, Ssol])
    NSsol = NSsol.reshape((-1,))

    return NSsol


# Jacobian (nNb+nNw+nS, nNb+nNw+nS) used by fsolve to calculate the steady-state (LI) characteristic
def Jac_cw_1D_transp(y, *args):

    SCHtransp, nNb, nNw, nS, Vr, c_act, c_injb, c_injw, It, c_diff, gln, c_st, Ntr, epsilon, GamZ, beta, tauNb, tauNw, tauCap, tauEsc, taub, tauw, tauS = args

    Nw = y[1:nNb+nNw]
    Nw = Nw[:, np.newaxis]
    S = y[nNb+nNw:]
    S = S[:, np.newaxis]

    Na = max(np.matmul(c_act, Nw), 1)                       # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # instantaneous linear gain coefficient, fitted to logarithmic gain function (-> time domain)
    #print(f'Na = {Na}; g0 = {g0}')

    diag_nNw = np.zeros((nNw, nNw), int)
    np.fill_diagonal(diag_nNw, 1)
    diag_nS = np.zeros((nS, nS), int)
    np.fill_diagonal(diag_nS, 1)
    slice_nSnNw = np.zeros((nS, nNw))
    slice_nSnNw[:, 0] = 1                                   # spontaneous emission proportional to Nw0 only (first column)

    Ngain = g0 * (np.squeeze(np.matmul(c_st, Nw)) - c_st[:, :, 0] * Ntr)
    Scompr = 1 + epsilon * S
    Sterm = S / Scompr
    #Stot = np.sum(S)
    #Sratio = S / Stot

    d00 = np.atleast_1d(-1/taub)[:, np.newaxis]                                             # dfNb/dNb (1,  1)
    d10 = c_injw / Vr / tauCap                                                              # dfNw/dNb (nNw,1)
    d20 = np.zeros((nS, 1))                                                                 # dfS/dNb  (nS, 1)
    
    d01 = np.zeros((1, nNw))                                                                # dfNb/dNw (1,  nNw)
    d01[0, 0] = Vr / tauEsc
    d11 = -(1/tauw + c_diff) * diag_nNw - g0 * np.sum(c_st * Sterm[:,:,np.newaxis] , 0)     # dfNw/dNw (nNw,nNw)
    d21 = GamZ * beta / tauNw * slice_nSnNw + GamZ * g0 * c_st[:, 0, :] * Sterm             # dfS/dNw  (nS, nNw)
    #21 = GamZ * beta / tauNw * slice_nSnNw * Sratio + GamZ * g0 * c_st[:, 0, :] * Sterm    # dfS/dNw  (nS, nNw)
    
    d02 = np.zeros((1, nS))                                                                 # dfNb/dS  (1,  nS)
    d12 = -(Ngain / Scompr**2).T                                                            # dfNw/dS  (nNw,nS)
    d22 = (-1 / tauS + GamZ * Ngain[:, 0] / Scompr**2) * diag_nS                            # dfS/dS   (nS, nS)
    #d22 = (-1 / tauS + GamZ * beta * Nw[0, :] / tauNw / Stot**2 * (diag_nS * Stot - np.tile(S.T, (nS, 1))) + GamZ * Ngain[:, 0] / Scompr**2) * diag_nS

    col0 = np.concatenate((d00, d10, d20), 0)
    col1 = np.concatenate((d01, d11, d21), 0)
    col2 = np.concatenate((d02, d12, d22), 0)

    return np.concatenate((col0, col1, col2), 1)


# system of ODEs to be solved for dynamic response calculation using solve_ivp solver
def ODEsolver_1D_transp(t, y, *args):

    SCHtransp, nNb, nNw, nS, Vr, c_act, c_injb, c_injw, It, c_diff, gln, c_st, Ntr, epsilon, GamZ, beta, tauNb, tauNw, tauCap, tauEsc, taub, tauw, tauS = args

    Nb = y[0, 0]
    Nw = y[1:nNb+nNw, 0]
    Nw = Nw[:, np.newaxis]
    S = y[nNb+nNw:, 0]
    S = S[:, np.newaxis]

    Na = max(np.matmul(c_act, Nw), 1)                       # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # instantaneous linear gain coefficient, fitted to logarithmic gain function (-> time domain)

    Injb = c_injb * It(t)                                   # current injection
    
    if SCHtransp == True:
        Cap = Nb / tauCap                                   # carrier capture term from the barriers into the QWs
        Esc = Nw / tauEsc                                   # carrier escape term from the QWs into the barriers
        Rnstb = Nb / tauNb                                  # non-stimulated carrier recombination in the barriers
    else:  
        Cap = Injb                                          # carrier capture term from the barriers into the QWs
        Esc = np.zeros((nNw, 1))                            # carrier escape term from the QWs into the barriers
        Rnstb = 0                                           # non-stimulated carrier recombination in the barriers

    Injw = c_injw / Vr * Cap                                # current injection into the barriers
    Rnstw = Nw / tauNw                                      # non-stimulated carrier recombination
    Diff = c_diff * Nw                                      # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, Nw)) \
          - c_st[:, :, 0] * Ntr) * S / (1 + epsilon * S)    # stimulated recombination
    Rolo = S / tauS                                         # optical losses and outcoupling term

    dNbdt = Injb - Rnstb - Cap + Vr * Esc[0,:][:, np.newaxis]
    dNwdt = Injw - Rnstw - Diff - Esc - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rolo + GamZ * beta * Rnstw[0, :] + GamZ * Rst[:, 0][:, np.newaxis]
    
    return np.concatenate([dNbdt, dNwdt, dSdt])


# Jacobian (nNb+nNw+nS, nNb+nNw+nS) used by solve_ivp to calculate the dynamic response
def Jac_1D_transp(t, y, *args):

    SCHtransp, nNb, nNw, nS, Vr, c_act, c_injb, c_injw, It, c_diff, gln, c_st, Ntr, epsilon, GamZ, beta, tauNb, tauNw, tauCap, tauEsc, taub, tauw, tauS = args
    

    Nw = y[1:nNb+nNw]
    Nw = Nw[:, np.newaxis]
    S = y[nNb+nNw:]
    S = S[:, np.newaxis]

    Na = max(np.matmul(c_act, Nw), 1)                       # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # instantaneous linear gain coefficient, fitted to logarithmic gain function (-> time domain)
    #print(f'Na = {Na}; g0 = {g0}')

    diag_nNw = np.zeros((nNw, nNw), int)
    np.fill_diagonal(diag_nNw, 1)
    diag_nS = np.zeros((nS, nS), int)
    np.fill_diagonal(diag_nS, 1)
    slice_nNwnS = np.zeros((nS, nNw))
    slice_nNwnS[:, 0] = 1                                   # spontaneous emission proportional to Nw0 only (first column)

    Ngain = g0 * (np.squeeze(np.matmul(c_st, Nw)) - c_st[:, :, 0] * Ntr)
    Scompr = 1 + epsilon * S
    Sterm = S / Scompr

    d00 = np.atleast_1d(-1/taub)[:, np.newaxis]                                             # dfNb/dNb (1,  1)
    d10 = c_injw / Vr / tauCap                                                              # dfNw/dNb (nNw,1)
    d20 = np.zeros((nS, 1))                                                                 # dfS/dNb  (nS, 1)
    
    d01 = np.zeros((1, nNw))                                                                # dfNb/dNw (1,  nNw)
    d01[0, 0] = Vr / tauEsc
    d11 = -(1/tauw + c_diff) * diag_nNw - g0 * np.sum(c_st * Sterm[:,:,np.newaxis] , 0)     # dfNw/dNw (nNw,nNw)
    d21 = GamZ * beta / tauNw * slice_nNwnS + GamZ * g0 * c_st[:, 0, :] * Sterm             # dfS/dNw  (nS, nNw)
    
    d02 = np.zeros((1, nS))                                                                 # dfNb/dS  (1,  nS)
    d12 = -(Ngain / Scompr**2).T                                                            # dfNw/dS  (nNw,nS)
    d22 = (-1 / tauS + GamZ * Ngain[:, 0] / Scompr**2) * diag_nS                            # dfS/dS   (nS, nS)

    col0 = np.concatenate((d00, d10, d20), 0)
    col1 = np.concatenate((d01, d11, d21), 0)
    col2 = np.concatenate((d02, d12, d22), 0)

    return np.concatenate((col0, col1, col2), 1)


# system of ODEs to be solved for dynamic response calculation using finite differences
def FDsolver_1D_transp(SCHtransp, Nbto, Nwto, Sto, nNw, Vr, c_act, c_injb, c_injw, It, c_diff, gln, c_st, Ntr, epsilon, GamZ, beta, tauNb, tauNw, tauCap, tauEsc, tauS):
    
    Na = max(np.matmul(c_act, Nwto), 1)                     # average carrier density over the active area
    g0 = gln * np.log(Na / Ntr) / (Na - Ntr)                # instantaneous linear gain coefficient, fitted to logarithmic gain function (-> time domain)

    Injb = c_injb * It                                      # current injection into the barriers

    if SCHtransp == True:
        Cap = Nbto / tauCap                                 # carrier capture term from the barriers into the QWs   (1,  1)
        Esc = Nwto / tauEsc                                 # carrier escape term from the QWs into the barriers    (nNw,1)
        Rnstb = Nbto / tauNb                                # non-stimulated carrier recombination in the barriers  (1,  1)
    else:  
        Cap = Injb                                          # carrier capture term from the barriers into the QWs   (1,  1)
        Esc = np.zeros((nNw, 1))                            # carrier escape term from the QWs into the barriers    (nNw,1)
        Rnstb = 0                                           # non-stimulated carrier recombination in the barriers  (1,  1)
        
    Injw = c_injw / Vr * Cap                                # carrier injection into the QWs
    Rnstw = Nwto / tauNw                                    # non-stimulated carrier recombination in the QWs
    Diff = c_diff * Nwto                                    # carrier diffusion
    Rst = g0 * (np.squeeze(np.matmul(c_st, Nwto)) \
        - c_st[:, :, 0] * Ntr) * Sto / (1 + epsilon * Sto)  # stimulated recombination
    Rsp = GamZ * beta * Rnstw[0, :]                         # spontaneous emission term
    Rolo = Sto / tauS                                       # optical losses and outcoupling term

    dNbdt = Injb - Rnstb - Cap + Vr * Esc[0,:]
    dNwdt = Injw - Rnstw - Diff - Esc - np.sum(Rst, 0)[:, np.newaxis]
    dSdt = -Rolo + Rsp + GamZ * Rst[:, 0][:, np.newaxis]
    
    return dNbdt, dNwdt, dSdt, Rsp

    
def HfResp(S, n, dt, fmaxplot, Parasitics, Hp):

    # 0. frequency vector (GHz)
    #f2 = np.arange(start=0, stop=fmaxplot, step=1/dt) * 1e-9    # cap f array at fmaxplot (default: 25GHz)
    
    # 0. Prep
    df = 1 / n / dt     # frequency resolution
    fmax = df * n / 2
    f = np.arange(start=0, stop=fmax, step=df)
    #f2 = fftfreq(n, dt)
    f = f[f<=fmaxplot] * 1e-9  # GHz
    # 1. response to small signal step
    P = np.sum(S, 0)
    # 2. derivative of unit step = impulse: H(f) = Y(f) = FFT(yimp(t))) = FFT(d/dt ystep(T)
    dP = np.gradient(P, dt)    # dP = np.diff(P) / dt as alternative
    # 3. FFT to transform to frequency domain (-> complex)
    DP = fft(dP)
    # 4. compute amplitude normalize, limite size to fmaxplot, convert to dB
    #H = np.sqrt(DP * np.conj(DP))
    H = np.abs(DP)
    H = H / H[0]
    H = H[:f.shape[0]]                  # cap H array at fmaxplot
    H = 10 * np.log10(H)
    # 6. similar steps for electrical parasitics response
    if Parasitics == True:
        Hp = np.abs(Hp)
        Hp = Hp[n//2:n//2+f.shape[0]]   # cap H array at fmaxplot
        Hp = 10 * np.log10(Hp)
    # (7. optional: filte0r to remove artefacts, e.g. in case rtol is too large in solve_ivp)
    #H = savgol_filter(H, window_length=(n/2), polyorder=4, mode='mirror')
    
    return f, H, Hp


def RINcalc(P, n, dt, nSeg, fmaxplot):

    # 0. frequency vector from 0 to fmaxplot (GHz)
    n = n // nSeg # time-domain response sliced into nSeg segments of length n/nSeg
    df = 1 / n / dt     # frequency resolution
    fmax = df * n / 2
    f = np.arange(start=0, stop=fmax, step=df)
    f = f[f<=fmaxplot] * 1e-9  # GHz
    #f = np.arange(start=0, stop=fmaxplot, step=1/n/dt)*1e-9    # cap f array at fmaxplot (default: 25GHz)
    # 1. power vector
    Pmean =  np.mean(P) # average over the complete vector (all segments)
    P2D = np.reshape(P, (nSeg, int(P.shape[0]/nSeg)))   # slice vector into nSeg segments
    # 2. FFT to transform to frequency domain (-> complex)
    DP = fft(P2D - Pmean, axis=1)
    # 3. compute amplitude as product of DP and its complex conjugate, remove 0 complex element
    DP = DP * np.conj(DP) / n
    DP = DP.real
    # 4. compute average over nSeg segments
    DP = np.mean(DP, 0)
    # 5. compute RIN in dB   
    RIN = 10 * np.log10(dt * DP / Pmean**2)

    return f[1:], RIN[1:f.shape[0]]


def parasiticsFilter(dt, n, Cp, Rm, Ca, Ra, It):
    fmax = 1 / dt / 2                           # 100GHz / 2 = 50GHz for dt = 10ps
    f = np.linspace(-fmax, fmax, n)             # symetrical frequency spectrum (full FFT)
    w = 2 * np.pi * f                           # angular frequency
    Hp = (1 / Cp / Rm / Ca / Ra) / (-w**2 - 1j * w * (Ca * Ra + Cp * Rm + Cp * Ra) / Cp / Ca / Ra / Ra + 1/ Cp / Ca / Rm / Ra)  # frequency response of passive parasitics circuit
    IHt = np.abs(ifft(Hp * fftshift(fft(It))))  # 1. It spectrum, 2. multiply by filter, 3. convert back to time domain, 4. take the absolute value
    # plt.plot(It)
    # plt.plot(IHt)
    # plt.show()

    return np.abs(Hp), IHt  # returns the real amplitude of Hp so that it can be saved in a json file


def VISTAS1D(sp, vp):
    print()

    # 1. key parameters --------------------------------------------------------------------------------
    # physical constants
    q = 1.602e-19               # elementary charge [C]
    h = 6.626e-34               # Planck constant [Js]
    #kB = 1.381e-23             # Boltzmann constant [J/K]
    c0 = 2.998e8                # vacuum speed of light [m/s]

    hvl = h * c0 * 1e9          # photon energy [J*m]
    vg = 100 * c0 / vp['ng']    # group velocity [cm/s]

    # cavity geometry parameters
    rCav = 3 * vp['rox']                        # optical field and carriers extend beyond the active layer, delimited by the oxide
    Vw = np.pi * rCav**2 * vp['nw'] * vp['dw']  # cavity volume (cm-3)
    Vb = np.pi * rCav**2 * vp['db']             # barriers volume (cm-3)
    Vr = Vw / Vb                                # ratio of QWs to barriers volume
    GamZ = vp['dw'] * vp['nw'] / vp['Leff']     # longitudinal optical confinement factor 
    nrho = 100  # radial steps
    nphi = 200  # azimuthal steps
    rho = np.linspace(0, rCav, nrho)
    rho = rho[:, np.newaxis]                    # np rank one array -> (1, nrho) vector (rank 2 array)
    phi = np.linspace(0, 2*np.pi, nphi)
    phi = phi[:, np.newaxis]                    # np rank one array -> (1, nphi) vector (rank 2 array)
    
    # misc parameters
    nNb = 1             # number of equations describing Nb, constant
    nNw = sp['nNw']     # single variable to improve readability
    if sp['SCHtransp'] == True:
        taub = (1/vp['tauNb'] + 1/vp['tauCap'])**-1     # "equivalent aggregated" carrier lifetime in the barriers
        tauw = (1/vp['tauNw'] + 1/vp['tauEsc'])**-1     # "equivalent aggregated" carrier lifetime in the QWs
    else:
        taub = vp['tauNb']
        tauw = vp['tauNw']

    # 2. normalized modes intensity profiles -----------------------------------------------------------

    # field profiles
    tStart = time.time()    
    nS, lvec, LPlm, ur, wl = LP_modes(vp['wl0']* 1e-9, vp['nc'], vp['dn'], vp['rox'] * 1e-2, nrho, rho.T * 1e-2, alpha = 10)

    # normalized intensity profiles
    Ur = np.square(ur)  # field amplitude -> intensity
    norm = np.trapz(Ur * rho.T, rho.T) * 2 / rCav**2
    norm = norm[:, np.newaxis]
    Ur = Ur / norm      # normalized intensity profile

    # size correction of parameters defined individually for each mode
    alpham = 1 / vp['Leff'] * np.log(1 / np.sqrt(vp['Rt'] * vp['Rb'])) * np.ones((nS, 1))   # mirror losses
    alphai = vp['alphai'] * np.ones((nS, 1))
    epsilon = vp['epsilon'] * np.ones((nS, 1))
    beta = vp['beta'] / nS * np.ones((nS, 1))   # division by nS to spread the spontaneous emission homogeneously over all modes (-> not additive!)
    GamZ = GamZ * np.ones((nS, 1))
    tEnd = time.time()
    print(f'Mode profiles computation ({nS} modes): {np.round(tEnd - tStart, 3)}s')

    # 3. radial injection current profile --------------------------------------------------------------

    iact = (rho <= vp['rox']) * 1   # "active" area, typically delimited by the oxide aperture
    ispread = (rho > vp['rox']) * np.exp(-(rho - vp['rox']) / vp['rs'])
    irho = iact + ispread
    norm = np.trapz(irho * rho, rho, axis = 0) * 2 / rCav**2
    irho = irho / norm
    
    # 4. coefficient arrays ----------------------------------------------------------------------------

    tStart = time.time()     
    # carrier radial series expansion terms J0i = J0(gammai * rho / rCav)
    gammai = np.asarray(jn_zeros(1, nNw-1))     # scipy function that calculates the roots of the Bessel function (excl first trivial 0 root!)
    gammai = np.insert(gammai, 0, 0, axis = 0)  # prepend 0 root, not included by jn_zero
    gammai = gammai[:, np.newaxis]
    J0i = jn(0, gammai * rho.T / rCav)          # bessel expansion terms J0i (nNw of them)

    # coefficients array for calculating the average carrier density Na over the active area Na
    c_act = np.zeros((1, nNw))
    c_act = 2 / vp['rox']**2 * np.trapz(J0i[:, 0:round(nrho * vp['rox'] / rCav)].T * rho[0:round(nrho * vp['rox'] / rCav), :], rho[0:round(nrho * vp['rox'] / rCav), :], axis = 0)    # integrates in the active area only
    c_act = c_act[:, np.newaxis].T  # (1, nNw)

    # coefficients array that captures the spatial overlap between current and carrier profiles
    c_injw = np.zeros((nNw, 1))
    c_injw = np.trapz(J0i.T * irho * rho, rho, axis = 0)    # (nNw, )
    c_injw = c_injw[:, np.newaxis]                          # (nNw, 1)
    c_injw = c_injw * 2 / (rCav * jn(0, gammai))**2         # (nNw, 1)

    c_injb =  vp['etai'] / q / Vb

    # coefficients array that captures the overlap between carrier expansion terms to model diffusion
    c_diff = np.zeros((nNw, 1))
    c_diff = vp['DN'] * (gammai / rCav)**2     # (nNw, 1)

    # coefficients array that captures the overlap between carrier expansion terms and mode profiles
    c_st = np.zeros((nS, nNw, nNw))
    for m in range(nS):
        for i in range(nNw):
            for j in range(nNw):
                    prodTmp = J0i[j, :] * Ur[m, :] * J0i[i, :] * rho.reshape((-1,)) # reshape -> dim (nrho, )
                    prodTmp = prodTmp[:, np.newaxis]
                    c_st[m, i, j] = 2 * vp['GamR'] * vg / (rCav * jn(0, gammai[i]))**2 * np.trapz(prodTmp, rho, axis = 0)    
    
    tEnd = time.time()
    print(f'Overlap parameters computation: {np.round(tEnd - tStart, 3)}s') 

    # 5. temp optical parameters derivation ------------------------------------------------------------
    
    tauS = 1 / vg / (alpham + alphai) # photon lifetime (to calculate the Rolo - the optical losses and outcoupling)
    F = (1 - vp['Rt']) / ((1 - vp['Rt']) + np.sqrt(vp['Rt'] / vp['Rb']) * (1 - vp['Rb'])) # fraction of power emitted at the top facet
    etaOpt = F * alpham / (alphai + alpham)             # optical efficiency (etaOpt*etai=etaD)  
    S2P = etaOpt * hvl / vp['wl0'] * Vw / tauS / GamZ * 1e3 # conversion factor photon density -> optical power (nS, 1)
    
    # 6. current and time characteristics --------------------------------------------------------------

    tStart = time.time()
    Hp = np.zeros((1,1))                # initialization as empty np array to allow saving results even if no parasitics filtering is calculated
    if sp['odeSolver'] == 'Finite Diff.':
        dt = sp['dtFD']
    else:
        dt = sp['dt']

    teval = np.arange(0, sp['tmax'], dt)
    n = len(teval)

    minGrad = (sp['Ion']-sp['Ioff']) * dt * 1e7 # min change of signal defining beginning of pulse (for step and pulse), accounting for different timesteps

    if sp['modFormat'] == 'step':
        It = np.zeros((2 * n))  # current vector of twice the length (FFT requires symetrical signal around zero)
        It[n:] = sp['Ion']
        if sp['Parasitics'] == True:
            Hp, It = parasiticsFilter(dt, 2*n, vp['Cp'], vp['Rm'], vp['Ca'], vp['Ra'], It)
        
        # quick&dirty fix for unexplained time-domain shift of the filtered signal (FFT phase?)
        sStart = np.where(np.gradient(It) > minGrad)    # sStart is a tuple -> sStart[0] returns an array, sStart[0][0] the first element of that array
        It = It[sStart[0][0]:sStart[0][0] + n]
        It[0] = 0                                       # to calculate NSinit@I(t0)=0

    elif sp['modFormat'] == 'pulse':
        It = np.ones((2 * n)) * sp['Ioff']              # current vector of twice the length (FFT requires symetrical signal around zero)
        It[n//2:n] =  sp['Ion']
        if sp['Parasitics'] == True:
            Hp, It = parasiticsFilter(dt, 2*n, vp['Cp'], vp['Rm'], vp['Ca'], vp['Ra'], It)
        
        # quick&dirty fix for unexplained time-domain shift of the filtered signal (FFT phase?)
        sStart = np.where(np.gradient(It)>minGrad)      # sStart is a tuple -> sStart[0] returns an array, sStart[0][0] the first element of that array
        It = It[sStart[0][0]:sStart[0][0] + n]
        It[0] = sp['Ioff']                              # to compute NSinit@I(t0)=Ioff

    elif sp['modFormat'] == 'random bits':
        nb = int(sp['tmax'] / sp['tb']) # number of bits
        sequence = np.random.randint(2, size = nb)
        cb = int(sp['tb'] / dt)         # number of time steps for one bit
        It = np.ones((n)) * sp['Ion']   # current vector
        for i in range(nb):
            It[i * cb : (i + 1) * cb] = It[i * cb : (i + 1) * cb] - (sp['Ion'] - sp['Ioff']) * sequence[i]
        if sp['Parasitics'] == True:
            Hp, It = parasiticsFilter(dt, n, vp['Cp'], vp['Rm'], vp['Ca'], vp['Ra'], It)

    elif sp['modFormat'] == 'small signal':
        if sp['Hfplot'] == True:        # tmax disabled and calculated automatically
            n = sp['nH'] + 1            # e.g. nH = 2**14 = 16'384 points -> df = 1/dt/n = 1/tmax = 62 MHz for dt = 1ps (typical for FD) or 6mHz for dt = 10ps
            sp['tmax'] = n * dt
        teval = np.arange(0, sp['tmax'], dt)
        It = np.ones((n)) * (sp['Ion'] + sp['Iss'])
        It[0] = sp['Ion']               # to compute NSinit@I(t0)=Ion
        if sp['Parasitics'] == True:
            Hp, _ = parasiticsFilter(dt, n, vp['Cp'], vp['Rm'], vp['Ca'], vp['Ra'], It) # in the small signal case, the signal is not affected to depict the intrinsic frequency response (parasitics contribution added separately on the plot)

    elif sp['modFormat'] == 'steady state':
        if sp['RINplot'] == True:       # tmax disabled and calculated automatically
            n = sp['nSeg'] * sp['nRIN'] # 8 segments of length e.g. nRIN = 2**15 = 32'768 each  -> 8* 32.75ns = 262ns, df = 1/dt/n = 1/tmax = 31 MHz (for dt = 1ps)
            sp['tmax'] = n * dt
        teval = np.arange(0, sp['tmax'], dt)
        It = np.ones((n)) * sp['Ion']   # current vector
    tEnd = time.time()
    print(f'Current signal generation: {np.round(tEnd - tStart, 3)}s')

       
    # 7. steady-state (LI) solution --------------------------------------------------------------------

    tStart = time.time() 
    Imax = 10e-3                                    # current range for LI characteristic
    dI = 0.01e-3                                    # current steps (must be small enough to ensure convergence)

    Icw = np.arange(0, Imax, dI)                    # current vector
    NScw = np.zeros((nNb+nNw+nS, Icw.shape[0]))     # NS solution matrix (nNb+nNw+nS, Imax / dI)

    # analytical solution (below threshold) as "seed"
    A = np.zeros((nNb+nNw+nS, nNb+nNw+nS))          # matrix to solve analytically system of linear Nb, Nw and Sm equations below threshold
    b = np.zeros((nNb+nNw+nS, 1))                   # b vector (Ax = b)
    
    for i in range(nNw):
        A[i+1, i+1] = -1/tauw - c_diff[i, 0]
        A[nNb+nNw:, 1] = -(GamZ * beta / vp['tauNw']).reshape((-1,))
    for i in range(nS):
        A[nNb+nNw+i, nNb+nNw+i] = 1 / tauS[i]

    if sp['SCHtransp'] == True:
        A[0, 0] = 1/taub                                                    # coefficients for Nb equation
        A[1:nNw+1, 0] = c_injw.reshape((-1,)) / Vr / vp['tauCap']           # coefficients for Nw equations
        A[0, 1] = -Vr / vp['tauEsc']
        b[0, 0] = c_injb * Icw[1]
        NScw[:, 1] = np.linalg.solve(A, b).reshape((-1,))                   # analytical solution for 2nd current step, calculated through matrix inversion
    else:
        b[1:nNb+nNw,0] = -c_injw.reshape((-1,)) * c_injb / Vr * Icw[1]
        NScw[1:, 1] = np.linalg.solve(A[1:, 1:], b[1:,:]).reshape((-1,))    # first row and column (related to Nb) removed from A and b matrices

    # main solver loop: sol@I0: 0, sol@I1: NScw[:, 1] calculated analytically above, sol@I3...n calculated below using fsolve and Jacobian
    for i in range(2, Icw.shape[0], 1):
        args = (sp['SCHtransp'], nNb, nNw, nS, Vr, c_act, c_injb, c_injw, Icw[i], c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, GamZ, beta, vp['tauNb'], vp['tauNw'], vp['tauCap'], vp['tauEsc'], taub, tauw, tauS)
        NScw[:, i] = fsolve(cw_1D_transp, NScw[:, i-1], args = args, fprime = Jac_cw_1D_transp)

    # interpolation function for later use
    NScwInterp = interp1d(x = Icw, y = NScw, fill_value = "extrapolate")    # extrapolates the LI characteristic -> NScwInterp(I)
    tEnd = time.time()
    print(f'LI solution: {np.round(tEnd - tStart, 3)}s')

    # 8. solver main loop ------------------------------------------------------------------------------

    tStart = time.time() 
    Nb = np.zeros((nNb, n))            # Nb(t)
    Nw = np.zeros((nNw, n))            # Nw0(t), Nw1(t), ..., NnNw(t)
    S = np.zeros((nS, n))              # multimode photon number matrix

    FS = np.zeros((nS, 1))              # modal noise terms

    NSinit = NScwInterp(It[0])          # result of steady-state (LI) characteristic (extrapolated) as starting point to avoid initial fluctuations
    
    f = np.zeros((1,1))                 # initialization as empty np array to allow saving results even if no freq response or RIN is calculated
    H = np.zeros((1,1))                 # initialization as empty np array to allow saving results even if no freq response is calculated
    RIN = np.zeros((1,1))               # initialization as empty np array to allow saving results even if no RIN is calculated

    if sp['odeSolver'] == 'Finite Diff.':   # 8a. solution of system of ODEs using finite differences with constant time-step

        Nbto = np.zeros((nNb, 1))       # temp variable passed to FD calculation function ('o' for 'old')
        Nwto = np.zeros((nNw, 1))       # temp variable passed to FD calculation function ('o' for 'old')
        Sto = np.zeros((nS, 1))         # temp variable passed to FD calculation function ('o' for 'old')

        Nb[:, 0] = NSinit[0]            # first point initialized at the LI value: Nb @ I(t=0)
        Nw[:, 0] = NSinit[1:nNb+nNw]    # first point initialized at the LI value: Nw @ I(t=0)
        S[:, 0] = NSinit[nNb+nNw:]      # first point initialized at the LI value: Sm @ I(t=0)
     
        for ti in range(n - 1):
            Nbto[0, 0] = Nb[0, ti]
            Nwto[:, 0] = Nw[:, ti]
            Sto[:, 0] = S[:, ti]

            dNbdt, dNwdt, dSdt, Rsp = FDsolver_1D_transp(sp['SCHtransp'], Nbto, Nwto, Sto, nNw, Vr, c_act, c_injb, c_injw, It[ti], c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, GamZ, beta, vp['tauNb'], vp['tauNw'], vp['tauCap'], vp['tauEsc'], tauS) # rhs of system of ODEs
            
            Nbtn = Nbto + dt * dNbdt            # Finite Differences step ('n' in Nbtn stands for for 'new')
            Nwtn = Nwto + dt * dNwdt            # Finite Differences step          
            Stn = Sto + dt * dSdt               # Finite Differences step

            # addition of noise
            if sp['Noise'] == True:
                xS = np.random.randn(nS, 1)                                     # random variable with zero mean and unit stdev
                FS = np.sqrt(np.maximum(2 * Rsp * Sto / dt, 0))* xS             # modal noise terms
                Stn =np.maximum(Stn + dt * FS, 0)                               # updated photon density
                
                xNw = np.random.randn(1)                                                            # random variable with zero mean and unit stdev
                FNw= -np.sum(FS) + np.sqrt(2 * np.maximum(Nwto[0, 0], 0) / vp['tauNw'] / dt) * xNw  # Nw00 noise term
                Nwtn[0, 0] = np.maximum(Nwtn[0, 0] + dt * FNw, 0)                                   # updated carrier density Nw00

                xNb = np.random.randn(1) * int(sp['SCHtransp'] == True)         # random variable with zero mean and unit stdev (forced to 0 if transport not simulated)
                FNb=np.sqrt(2 * np.maximum(Nbto, 0) / vp['tauNb'] / dt) * xNb   # Nw00 noise term
                Nbtn = np.maximum(Nbtn + dt * FNb, 0)                           # updated carrier density Nw00
            
            Nb[0, ti+1] = Nbtn.reshape((-1,))
            Nw[:, ti+1] = Nwtn.reshape((-1,))
            S[:, ti+1] = Stn.reshape((-1,))

        tEnd = time.time()
        print(f'FD ODE solution: {np.round(tEnd - tStart, 3)}s')

        if sp['Hfplot'] == True:   # must be calculated prior to subsampling
            tStart = time.time()
            f, H, Hp = HfResp(S[:, 1:], n-1, dt, sp['fmaxplot'], sp['Parasitics'], Hp)  # compute frequency response
            tEnd = time.time()
            print(f'Frequency response calculation: {np.round(tEnd - tStart, 3)}s')

        if sp['RINplot'] == True:
            tStart = time.time()
            f, RIN = RINcalc(np.sum(S* S2P, 0), n, dt, sp['nSeg'], sp['fmaxplot'])      # compute RIN spectrum
            tEnd = time.time()
            print(f'RIN calculation: {np.round(tEnd - tStart, 3)}s')

        # small timestep dtFD needed to ensure convergence, but too dense for post-processing -> subsampling to larger timestep dt  
        Nb = Nb[:, ::int(sp['dt']/dt)]
        Nw = Nw[:, ::int(sp['dt']/dt)]
        S = S[:, ::int(sp['dt']/dt)]
        teval = teval[::int(sp['dt']/dt)]
        It = It[::int(sp['dt']/dt)]
        n = len(teval)
        dt = sp['dt']
 
    else:   # 8b.solution of system of ODEs using 'solve_ivp'

        tStart = time.time()        
        Itinterp = interp1d(x = teval, y = It, fill_value = "extrapolate") # current changes over time and must match the solve_ivp integration points
        
        args = (sp['SCHtransp'], nNb, nNw, nS, Vr, c_act, c_injb, c_injw, Itinterp, c_diff, vp['gln'], c_st, vp['Ntr'], epsilon, GamZ, beta, vp['tauNb'], vp['tauNw'], vp['tauCap'], vp['tauEsc'], taub, tauw, tauS)       
        sol = solve_ivp(ODEsolver_1D_transp, (0, sp['tmax']), NSinit, t_eval=teval, method=sp['odeSolver'], dense_output=True, vectorized=True, args=args, rtol=1e-6, atol=1e-6, jac=Jac_1D_transp)

        Nb = sol.y[0]
        Nw = sol.y[1:nNb+nNw]
        S = sol.y[nNb+nNw:]

        tEnd = time.time()
        print(f'solve_ivp ODE solution: {np.round(tEnd - tStart, 3)}s')

        if sp['modFormat'] == 'small signal':
            tStart = time.time()
            f, H, Hp = HfResp(S[:, 1:], n-1, dt, sp['fmaxplot'], sp['Parasitics'], Hp)    # compute frequency response
            tEnd = time.time()
            print(f'Frequency response calculation: {np.round(tEnd - tStart, 3)}s')

    # 9. output dictionary sr constructed from individual variables 
    
    sr = {
        "phi": phi,
        "nphi": nphi,
        "rho": rho,
        "nrho": nrho,
        "nNw": nNw,
        "J0i": J0i,
        "nS": nS,
        "LPlm": LPlm,
        "lvec": lvec,
        "Ur": Ur,
        "wl": wl,
        "Icw": Icw,
        "NScw": NScw,
        "It": It,
        "teval": teval,
        "f": f,
        "H": H,
        "Hp": Hp,
        "RIN": RIN,
        "S2P": S2P,
        "Nb": Nb,
        "Nw": Nw,
        "S": S,
        }

    return sr

    
def main():
    
    with open('last_params.json', 'r') as openfile:
        params = json.load(openfile)    # dictionary
    sp, vp = params['simParams'], params['vcselParams']
    del params

    sr = VISTAS1D(sp, vp)


if __name__ == "__main__":
    main()