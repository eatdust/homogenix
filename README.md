# Homogenix: parallel convolution by PSFex kernels

### Summary

Homogenix is a simple modern fortran code designed to perform
fast PSF homogenization from the
[psfex](https://github.com/astromatic/psfex) variable kernels.


### Compilation

Please ensure that you have a working installation of the **gfortran**
compiler (possibly its MPI wrappers, or alternatives) and the
[cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) library. Editing the
provided Makefile is certainly needed to specify the install location
of these libraries.

### Usage is Iraf inspired

Dealing with multiple files (MPI+OMP parallelized):

        homogenix @inlistfiles @inkernelfiles @outlistfiles

For a single image (OMP parallelized):

        homogenix inimage kernelcube outimage

---

