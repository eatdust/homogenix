!   This file is part of homogenix
!
!   Copyright (C) 2024 C. Ringeval
!   
!   gaialaxy is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   gaialaxy is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with homogenix.  If not, see <https://www.gnu.org/licenses/>.
!

module psfeed
  use precision,  only : fsp, fdp
  use iofits, only : read_kernel_homofits, read_key_homofits

  implicit none

  private
  
  real(fsp) :: order,xoffset, yoffset, xscale, yscale
  real(fsp), dimension(:,:,:), allocatable :: kernels
  integer :: na, nb, nk
  

  public initialize_kernels, check_kernels, free_kernels
  public get_psf_xsize, get_psf_ysize, get_psf_kernel
  
contains

 
 
  subroutine initialize_kernels(filename)
    implicit none
    character(len=*), intent(in) :: filename

    write(*,*)'reading kernels: ',filename
    call read_kernel_homofits(filename,kernels)
    na = size(kernels,1)
    nb = size(kernels,2)
    nk = size(kernels,3)
    
    if (nk.gt.6) stop 'initialize_kernels: only up to order 2 right now!'

    call read_key_homofits(filename,'POLZERO1',xoffset)
    call read_key_homofits(filename,'POLSCAL1',xscale)
    
    call read_key_homofits(filename,'POLZERO2',yoffset)
    call read_key_homofits(filename,'POLSCAL2',yscale)
            
    
  end subroutine initialize_kernels


  function get_psf_xsize()
    implicit none
    integer :: get_psf_xsize
    if (.not.check_kernels()) stop 'get_psf_xsize: kernel not found!'

    get_psf_xsize = na

  end function get_psf_xsize

  
  function get_psf_ysize()
    implicit none
    integer :: get_psf_ysize
    if (.not.check_kernels()) stop 'get_psf_ysize: kernel not found!'

    get_psf_ysize = nb

  end function get_psf_ysize

  
  function check_kernels()
    implicit none
    logical :: check_kernels

    check_kernels = allocated(kernels)

  end function check_kernels


  subroutine free_kernels()
    implicit none

    if (check_kernels()) then
       deallocate(kernels)
    endif

  end subroutine free_kernels
      

  
  subroutine get_psf_kernel(x,y,psf)
    implicit none
    real(fsp), intent(in) :: x,y
    real(fsp), dimension(:,:) :: psf
    
    real(fsp) :: xnorm, ynorm
    xnorm = (x-xoffset)/xscale
    ynorm = (y-yoffset)/yscale
    
    psf = kernels(:,:,1) + xnorm * kernels(:,:,2) + ynorm * kernels(:,:,3) &
         + xnorm*ynorm * kernels(:,:,4) + xnorm*xnorm * kernels(:,:,5) + ynorm*ynorm * kernels(:,:,6)
    
  end subroutine get_psf_kernel

  

end module psfeed
