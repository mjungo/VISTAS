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
from matplotlib import cm

import numpy as np


# plot 2D profile of each mode
def plotModes2D(Ur, LPlm, lvec, nS, rho, nrho, phi, nphi, nfig):
    # createmesh
    P, R = np.meshgrid(phi, rho)        # in polar coordinates
    X, Y = R * np.cos(P), R * np.sin(P) # in cartesian coordinates

    # format 2D profiles
    UR = np.tile(Ur, (nphi, 1, 1))      # radial modal intensities (nphi, nS, nrho)
    UR = np.transpose(UR, (1, 2, 0))    # radial modal intensities (nS, nrho, nphi)
    UPc, UPs = np.zeros((nS, nrho, nphi)), np.zeros((nS, nrho, nphi))   # azimuthal modal intensities ("c" for cos(l*phi) modes, "s" for sin(l*phi) modes)
    for m in range(nS):
        UPc[m, :, :], UPs[m, :, :] = np.cos(lvec[m] * P)**2, np.sin(lvec[m] * P)**2
    Uc, Us = UPc * UR, UPs * UR

    # plot normalized intensity profiles
    if nS == 1:
        fig = plt.figure(nfig)
        ax = plt.subplot(1, 1, 1, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[0, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[0] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        nfig = nfig + 1 # number of next figure

    elif nS == 2:
        fig = plt.figure(nfig)
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[0, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[0] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[1, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[1] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        nfig = nfig + 1 # number of next figure

    elif nS == 3:   # this case actually doesn't happen for LP modes (modes 3 and 4 with same cutoff frequency), but may be useful for different mode solvers
        fig = plt.figure(nfig)
        ax = fig.add_subplot(1, 3, 1, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[0, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[0] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(1, 3, 2, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[1, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[1] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(1, 3, 3, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[2, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[2] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        nfig = nfig + 1 # number of next figure

    elif nS == 4:
        fig = plt.figure(nfig)
        ax = fig.add_subplot(2, 2, 1, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[0, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[0] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(2, 2, 2, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[1, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[1] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(2, 2, 3, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[2, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[2] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])

        ax = fig.add_subplot(2, 2, 4, projection='3d')
        ax.plot_surface(X*1e6, Y*1e6, Uc[3, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
        ax.set_title(LPlm[3] + 'c' + '\nnormalized intensity (a.u.)')
        ax.set_xlabel('cavity x(\u03BCm)')
        ax.set_ylabel('cavity y(\u03BCm)')
        ax.zaxis.set_visible(False)
        ax.set_zticklabels([])    

        nfig = nfig + 1 # number of next figure

    else:   # generic case, nS > 4, plot on a 2x3 subplots matrix
        figs = nS // 6  # number of figures with 6 subplots each required to display all nS modes
        if nS % 6 != 0: figs = figs + 1 # add one "uncomplete" figure in case nS is not a multiple of 6
        for f in range(figs):
            fig = plt.figure(nfig)
            for m in range(6*f,min(6*f+6, nS), 1):  # on the last figure, stops at nS (uncomplete figure)
                ax = fig.add_subplot(2, 3, m-6*f+1, projection='3d')
                ax.plot_surface(X*1e6, Y*1e6, Uc[m, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
                ax.set_title(LPlm[m] + 'c')
                ax.set_xlabel('cavity x(\u03BCm)')
                ax.set_ylabel('cavity y(\u03BCm)')
                ax.zaxis.set_visible(False)
                ax.set_zticklabels([])
            nfig = nfig + 1 # number of next figure
    
    # for m in range(nS):   # this plots each mode profile in a separate figure
    #     ax = plt.subplot(projection='3d')
    #     ax.plot_surface(X*1e6, Y*1e6, Uc[m, :, :], antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
    #     #ax = plt.subplot()
    #     #ax.pcolormesh(X*1e6, Y*1e6, Uc[m, :, :], antialiased=False, cmap=cm.seismic)
    #     #plt.axis('scaled')
    #     #fig_title = LPlm[m] + 'c, ' + 'normalized intensity (a.u.)'
    #     ax.set_title(LPlm[m] + 'c' + '\nnormalized intensity (a.u.)')
    #     ax.set_xlabel('cavity x(\u03BCm)')
    #     ax.set_ylabel('cavity y(\u03BCm)')
    #     ax.zaxis.set_visible(False)
    #     ax.set_zticklabels([])

    return nfig


# plot the output optical power (incl. modal split)
def plotPower(x, P, modes, xlabel, nfig):
    plt.figure(nfig)
    plt.plot(x, P.T)
    plt.legend(modes)
    plt.plot(x, np.sum(P, 0), 'k--')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(xlabel)
    plt.ylabel('optical output power (mW)')
    plt.grid()

    return nfig + 1


# plot 2D carrier and optical (aggregated over all modes) profiles within the cavity
def plotNwS2D(Nw, J0i, S, Ur, rho, nrho, phi, nphi, nfig):
    # create mesh
    P, R = np.meshgrid(phi, rho)        # in polar coordinates
    X, Y = R * np.cos(P), R * np.sin(P) # in cartesian coordinates

    # reconstruct and plot 2D carrier profile
    Nw2Dprof = np.tile(J0i, (nphi, 1, 1))    # carrier expansion terms (nphi, nNw, nrho)
    Nw2Dprof = np.transpose(Nw2Dprof, (1, 2, 0))  # carrier expansion terms (nNw, nrho, nphi)
    Nw2Dprof = np.sum(Nw[:, np.newaxis, np.newaxis] * Nw2Dprof,0)
    fig = plt.figure(nfig)
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.plot_surface(X*1e6, Y*1e6, Nw2Dprof, antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
    ax.set_title('carrier density')
    ax.set_xlabel('cavity x(\u03BCm)')
    ax.set_ylabel('cavity y(\u03BCm)')
    ax.zaxis.set_visible(False)
    ax.set_zticklabels([])

    # reconstruct and plot 2D optical field profile
    S2Dprof = np.tile(Ur, (nphi, 1, 1))      # radial modal intensities (nphi, nS, nrho)
    S2Dprof = np.transpose(S2Dprof, (1, 2, 0))    # radial modal intensities (nS, nrho, nphi)
    S2Dprof = np.sum(S[:, np.newaxis, np.newaxis] * S2Dprof,0)
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_surface(X*1e6, Y*1e6, S2Dprof, antialiased=False, cmap=cm.seismic) # seismic, magma, cividis /// plot_surface, countour3D, plot_wireframe
    ax.set_title('photon density')
    ax.set_xlabel('cavity x(\u03BCm)')
    ax.set_ylabel('cavity y(\u03BCm)')
    ax.zaxis.set_visible(False)
    ax.set_zticklabels([])

    return nfig + 1


# plot frequency response spectrum
def plotH(f, H, Parasitics, Hp, fmin, fmax, ymin, ymax, nfig):
    plt.figure(nfig)
    plt.plot(f, H, lw=1)
    if Parasitics == True:
        plt.plot(f, Hp, lw=1)
        plt.plot(f, H + Hp, 'k', lw=1.5) #log(x*y) = log(x) + log(y)
        plt.legend(['intrinsic response', 'parasitics response', 'total response'])
    else:
        plt.legend(['intrinsic response'])
    plt.plot([fmin, fmax], [-3, -3], 'k-', lw=1,dashes=[2, 2])
    plt.xlim(xmin=fmin, xmax=fmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('frequency response (dB)')
    plt.grid()
    fmin = fmin + 1 - 1

    return nfig + 1


# plot spectra (frequency response and RIN)
def plotRIN(f, RIN, fmin, fmax, ymin, ymax, nfig):
    plt.figure(nfig)
    plt.plot(f, RIN)
    plt.xlim(xmin=fmin, xmax=fmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('Relative Intensity Noise (dB/Hz)')
    plt.grid()

    return nfig + 1


# plot simple eye diagram
def plotEye(teval, dt, tb, P, It, nfig):
    fig = plt.figure(nfig)

    n = teval.shape[0]
    if max(teval) > 10e-9:
        Itplot = It[n//4:n//4+int(10e-9/dt)]     # takes a 5ns sample starting at tmax / 4
        tplot = teval[n//4:n//4+int(10e-9/dt)]   # takes a 5ns sample starting at tmax / 4
    else:
        Itplot = It
        tplot = teval
    ax = fig.add_subplot(2, 1, 1)
    ax.plot(tplot*1e9, Itplot*1e3)
    ax.set_title('10ns drive signal sample')
    ax.set_ylabel('drive current (mA)')

    P = np.sum(P, 0)    # sum over all modes -> Ptot
    tmax = max(teval)
    ctb = int(tb / dt)   # number of time steps for one "bit"
    nb = int(tmax / tb)
    tplot = np.linspace(0, 2*tb, 2*ctb)
    #tplot = np.arange(0, 2*tb, dt)
    ax = fig.add_subplot(2, 1,2)
    for i in range(1, nb-1):
        ax.plot(tplot*1e9, P[int((i-0.5)*ctb):int((i+1.5)*ctb)], 'b')  # displays 2 bit lengths, from -0.5 to 1.5 * tb
    ax.set_title(str(nb) + ' periods')
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('optical output power (mW)')

    return nfig + 1