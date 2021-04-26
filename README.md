# VISTAS Python

VISTAS is a 2D, multimode, time-domain model aimed at simulating the static and dynamic behaviour of Vertical Cavity Surface Emitting Lasers (VCSEL).

VISTAS accounts dynamically for the spatial interactions between the optical field and carrier distributions in the active layer. Effects such as spatial hole burning, modal gain competition, or carrier diffusion losses can therefore be studied in details.

The original VISTAS algorithm was published in the early 2000s [1][2]. It is based on 2D spatially-dependent rate equations, mathematically transformed so as to remove any explicit spatial dependency. The resulting model consists therefore in a system of Ordinary Differential Equations (ODEs). It can be implemented in a fully vectorized fashion and exhibits several orders of magnitude improvement in terms of computational efficiency compared to the original spatially-dependent rate equations.

The model can also easily be extended with a series of important effects affecting the performance of VCSEL-based optical links, such as carrier transport into the quantum wells, noise, thermal effects, etc.

The aim of the project hosted in this repository is to develop a modern implementation of VISTAS in Python. Contributions are warmly welcomed!


References:

[1] M. Jungo, "Spatiotemporal VCSEL Model for Advanced Simulations of Optical Links," in Series in Quantum Electronics, vol. 30, edited by H. Baltes, P. Günter, U. Keller, F. K. Kneubühl, W. Lukosz, H. Mechior, and M. W. Sigrist, 1st ed. Konstanz: Hartung-Gorre Verlag, 2003, ISBN 3-89649-831-2

[2] M. Jungo, D. Erni and W. Baechtold, "VISTAS: a comprehensive system-oriented spatiotemporal VCSEL model," ” IEEE J. Sel. Top. Quantum Electron., vol. 9, no. 3, pp. 939–948, May/Jun. 2003