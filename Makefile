#   This file is part of homogenix
#
#   Copyright (C) 2024 C. Ringeval
#   
#   homogenix is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   homogenix is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with homogenix.  If not, see <https://www.gnu.org/licenses/>.


ext=$(shell uname | cut -c1-3)

#FC=mpif90
FC=gfortran

FFLAGS = -O2 -mavx2 -fopenmp
#FFLAGS= -DMPISCHED -fopenmp-simd

INCLUDE= -I/usr/include/
LFLAGS= -L/usr/lib64 -lcfitsio


OBJS= precision.o iofits.o iohomo.o psfeed.o convolutor.o
OBJSMPI = fifo.o iofifo.o scheduler.o

homogenix.$(ext) : $(OBJS) $(OBJSMPI) homogenix.o
	$(FC) $(FFLAGS) $(OBJS) $(OBJSMPI) homogenix.o $(LFLAGS) -o $@

%.o: %.F08
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.F03
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.f03
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $<
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean:
	rm *.$(ext) *.o *.mod


