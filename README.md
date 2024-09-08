# Homogenix: parallel convolution of fits image by psfex kernels

### Summary

Homogenix is a very simple modern fortran code designed to perform
fastPSF homogeneisation from the
[PSFex](https://github.com/astromatic/psfex) variables kernels.


### Compilation

Please ensure that you have a working installation of the **gfortran** and
**gcc** compilers (or alternatives), the
[cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) and
[wcslib](https://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/)
libraries. Editing the provided Makefile is certainly needed to
specify the install location of these libraries.

### Usage is Iraf inspired

Dealing with multiple files (MPI parallelised):

	homogenix @inlistfiles @inkernelfiles @outlistfiles

For a single image (OMP parallelised):

        homogenix inimage kernelcube outimage

---

