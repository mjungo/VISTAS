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

def params():

    # general physical parameters ---------------------------------------------

    wl0 = 858                           # emission wavelength @300K (nm)
    Rt = .997 * np.ones((20, 1))        # top mirror reflectivity
    Rb = .999 * np.ones((20, 1))        # top mirror reflectivity
    alpha_int = 30 * np.ones((20, 1))   # internal losses (cm-1)
    tau_N = 1.6e-9                      # carrier-lifetime (s)
    beta = 2e-4 * np.ones((20, 1))      # spontaneous recombination coeff.
    Gamma = 0.05 * np.ones((20, 1))     # optical confinement factor

    # carrier transport and diffusion parameters ------------------------------

    tau_esc = 400e-12                   # thermionic emission lifetime (s)
    tau_cap = 45e-12                    # ambipolar diffusion time (s)
    eta_i = 0.9                         # current injection efficiency
    rs = .5e-4                          # current spreading coefficient (cm)
    DN = 12                             # ambipolar diffusion coeff. (cm2/s)

    # gain parameters ---------------------------------------------------------

    Ntr = 1.8e18                        # transparency carrier density (cm-3)
    gln = 1300                          # logarithmic gain coefficient (cm-1)
    #g0 = 5e-16                         # linear gain coefficient (cm2)
    epsilon = 5e-17 * np.ones((20, 1))  # gain compression factor

    # cavity parameters ------------------------------------------------

    r_ox = 4.5e-4                       # oxide aperture radius (cm)
    r_cav = 2.5 * r_ox                  # cavity radius (cm)
    nc = 3.6                            # core equivalent refractive index
    ng = 4.2                            # group refractive index
    delta_n = .0008                     # equivalent fractional refractive index change dn
    L = 900e-7                          # effective cavity length (cm)
    dqw = 8e-7                          # single QW thickness (cm)
    nw = 3                              # number of quantum wells
    db = 40e-7                          # SCH thickness (cm)
    V_cav = np.pi * r_cav**2 * nw * dqw # cavity volume (cm-3)
    t = 0.7e-6                          # oxide thickness (cm)  

    # thermal parameters ------------------------------------------------------

    Ne = 2e-8                           # diode voltage parameter
    Rth = 3000                          # thermal resistance (K/W)
    c = 0.35                            # specific heat (J/g/C)
    density = 5.36                      # mass density (g/cm3)
    #Cth = density*c*pi*R^2*nw*dqw      # thermal capacitance (J/K)
    dwlT = .06                          # temperature coeff. of lFP (nm/K)
    gwl0 = 848                          # gain peak wavelength at T = 300K (nm)
    dgwlT = .27                         # temperature coeff. of lp (nm/K)
    g_width = 40                        # gain profile FWHM (nm)
    Tref = 250                          # reference temperature (K)
    Il0 = .0006                         # leakage parameter (A)
    a0 = -700                           # leakage parameter
    a1 = 5.4e-17                        # leakage parameter
    a2 = 2.4e-19                        # leakage parameter
    a3 = -3.4e21                        # leakage parameter
    b1 = .5e16                          # T-dependence coeff. of Ntr (cm-3/K)

    # parasitics --------------------------------------------------------------

    Rq = 50                             # current source resistance (ohm)					
    Cq = .5e-12                         # current source capacitance (F)
    Rw = .4                             # bond-wire resistance (ohm)
    Lw = 1e-9                           # bond-wire inductance (H)
    Cp = 0.5e-12                        # pad capacitance (F)
    Rs = 20                             # Bragg reflectors resistance (ohm)
    Ra = 30                             # cavity resistance (ohm)
    Ca = .5e-12                         # cavity capacitance (F)

    # external cavity parameters ----------------------------------

    d_ext = 40e-4                       # distance laser facet - fiber (cm)
    L_ext = 20                          # external cavity length (cm)
    n_ext = 1.46                        # external medium refractive index
    #R_ext = ((1-next)/(1+next))**2     # external power reflectance1
    rwg = 25e-4                         # fiber core radius (cm)

    # finite differences parameters -------------------------------------------

    nrho = 100                          # radial steps
    nphi = 200                          # azimuthal steps

    return \
        (wl0, Rt, Rb, alpha_int, beta, Gamma, \
            tau_N, tau_esc , tau_cap, eta_i, rs, DN, \
                Ntr, gln, epsilon, \
                    r_ox, r_cav, nc, ng, delta_n, L, V_cav, db, \
                        nrho, nphi)