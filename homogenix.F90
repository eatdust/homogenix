!   This file is part of homogenix
!
!   Copyright (C) 2024 C. Ringeval
!   
!   homogenix is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   homogenix is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with homogenix.  If not, see <https://www.gnu.org/licenses/>.
!
program homogenix
  use precision, only : fdp, fsp
#if defined MPI
  use mpi
#endif
  use scheduler, only : initialize_scheduler, free_scheduler,scheduled_size
  use scheduler, only : start_scheduling, irq_scheduler, stop_scheduling
  use scheduler, only : restore_scheduler, scheduler_save_queue, check_saved_queue
  use iohomo, only : fitsfile, get_arg_files , free_fitsfile
  use iofits, only : read_image_fits, copy_hdr_fits, overwrite_image_fits
  use psfeed, only : initialize_kernels, free_kernels
  use convolutor, only : pix_convolve

  implicit none

#if defined MPI
  integer :: mpiCode
#endif

  integer, save :: mpiRank = 0
  integer, save :: mpiSize = 1
!restore all processes from previous run
  logical, parameter :: cpRestart = .false.
!save queues each time one element is done
  logical, parameter :: cpSave = .true.
!if zero, assume same mpiSize between runs.
  integer, save :: mpiPrevSize = 0


  type(fitsfile), dimension(:), allocatable :: listfits
    
  integer :: nx,ny,i
  real(fsp), dimension(:,:), allocatable :: image, omage

  integer :: ifile, nfiles
  character(len=:), allocatable :: filekernels 
  character(len=:), allocatable :: fileimage 
  character(len=:), allocatable :: fileomage


  logical, parameter :: display = .true.
  

  call get_arg_files(listfits)
  nfiles = size(listfits(1)%files,1)
  
  
#ifdef MPI
  call MPI_INIT(mpiCode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiRank,mpiCode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSize,mpiCode)
#endif

  if (mpiPrevSize.eq.0) mpiPrevSize = mpiSize

  if (cpRestart) then
     call restore_scheduler(mpiPrevSize)
  else
     call initialize_scheduler(nfiles)
  endif

  if ((cpSave).and.(.not.cpRestart)) then
     if (check_saved_queue(mpiRank)) stop 'previous saved files present!'
  endif

  do

     call start_scheduling(ifile)

     if (display) then
        write(*,*)
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)'RANK= ',mpiRank,'   QSIZE= ',scheduled_size()
        write(*,*)'IMGFILE: ',listfits(1)%files(ifile+1)
        write(*,*)'KERFILE: ',listfits(2)%files(ifile+1)
        write(*,*)'OUTFILE: ',listfits(3)%files(ifile+1)
        write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)
     end if

     fileimage  = trim(listfits(1)%files(ifile+1))
     filekernels= trim(listfits(2)%files(ifile+1))
     fileomage  = trim(listfits(3)%files(ifile+1))


     call read_image_fits(fileimage,image)
     nx = size(image,1)
     ny = size(image,2)
     allocate(omage(nx,ny))

     call initialize_kernels(filekernels)

     call pix_convolve(image,omage)

     call free_kernels()

     call copy_hdr_fits(fileimage,fileomage)

     call overwrite_image_fits(fileomage,omage)

     deallocate(image,omage)
     deallocate(fileimage,filekernels,fileomage)

     call irq_scheduler()

     if (cpSave) then
        call scheduler_save_queue(mpiRank)
     endif

     if (stop_scheduling()) exit

  enddo

  call free_scheduler()

#ifdef MPI
  write(*,*)'process on barrier: mpiRank= ',mpiRank
  call MPI_BARRIER(MPI_COMM_WORLD,mpiCode)
  call MPI_FINALIZE(mpiCode)
#endif

  if (mpiRank.eq.0) then
     write(*,*)'all files done!'
  
     do i=1,size(listfits,1)
        call free_fitsfile(listfits(i))
     enddo
  endif

  

end program homogenix


