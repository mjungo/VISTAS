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


# plot the output optical power (incl. modal split)
def plotPower(x, P, modes, xlabel):
    plt.plot(x, P.T)
    plt.legend(modes)
    plt.plot(x, np.sum(P, 0), 'k--')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(xlabel)
    plt.ylabel('optical output power (mW)')
    plt.grid()
    plt.show()


# plot frequency response
def plotH(f, H):

    plt.plot(f, H)
    plt.xlim(xmin=0, xmax=10)
    plt.ylim(ymin=-10, ymax=15)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('frequency response (dB)')
    plt.grid()
    plt.show()


# plot simple eye diagram
def plotEye(teval, dt, tb, P):
    P = np.sum(P, 0)
    tmax = max(teval)
    ctb = int(tb / dt)   # number of time steps for one "bit"
    nb = int(tmax / tb)
    tplot = np.linspace(0, 2*tb, 2*ctb)
    #tplot = np.arange(0, 2*tb, dt)
    for i in range(1, nb-1):
        plt.plot(tplot, P[int((i-0.5)*ctb):int((i+1.5)*ctb)], 'b')  # displays 2 bit lengths, from -0.5 to 1.5 * tb
    plt.title(str(nb) + ' periods')
    plt.ylabel('optical output power (mW)')
    plt.xlabel('time (ns)')
    plt.show()