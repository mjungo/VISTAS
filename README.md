# VISTAS Python

VISTAS is a 2D, multimode, time-domain model aimed at simulating the static 
and dynamic behaviour of Vertical Cavity Surface Emitting Lasers (VCSEL).

VISTAS accounts dynamically for the spatial interactions between the optical field 
and carrier distributions in the active layer. Such effects as spatial hole 
burning, modal gain competition, or carrier diffusion losses can therefore 
be studied in details.

The original algorithm behind VISTAS  was published in the early 2000s [1][2].
It is based on 2D spatially-dependent rate equations, mathematically transformed
so as to remove any explicit spatial dependency. The resulting numerical 
efficiency is improved by several orders of magnitude.

The model can easily be extended with a series of important effects, 
such as carrier transport into the quantum wells, noise, thermal effects, etc.

The aim of the project stored in this repository is to develop a modern implementation 
of VISTAS in Python. Contributions are warmly welcomed!


References:
[1] M. Jungo, "Spatiotemporal VCSEL Model for Advanced Simulations of Optical 
Links," in Series in Quantum Electronics, vol. 30, edited by H. Baltes, P. Günter,
U. Keller, F. K. Kneubühl, W. Lukosz, H. Mechior, and M. W. Sigrist, 
1st ed. Konstanz: Hartung-Gorre Verlag, 2003, ISBN 3-89649-831-2
[2] M. Jungo, D. Erni and W. Baechtold, "VISTAS: a comprehensive system-oriented 
spatiotemporal VCSEL model," ” IEEE J. Sel. Top. Quantum Electron., 
vol. 9, no. 3, pp. 939–948, May/Jun. 2003