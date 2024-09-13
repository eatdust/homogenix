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

module psfeed
  use precision,  only : fsp, fdp
  use iofits, only : read_kernel_homofits, read_key_homofits

  implicit none

  private
  
  real(fsp) :: xoffset, yoffset, xscale, yscale
  real(fsp), dimension(:,:,:), allocatable :: kernels
  integer :: na, nb, nk
  
  logical, parameter :: display = .false.
  
  public initialize_kernels, check_kernels, free_kernels
  public get_psf_xsize, get_psf_ysize, get_norm_psf_kernel
  
contains

 
 
  subroutine initialize_kernels(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(fsp), dimension(:,:,:), allocatable :: cube
    integer :: k
    
    if (display) write(*,*)'reading kernels: ',filename

    call read_kernel_homofits(filename,cube)

    na = size(cube,1)
    nb = size(cube,2)
    nk = size(cube,3)

!for contiguous memory access within get_psf    
    allocate(kernels(nk,na,nb))
    do k=1,nk
       kernels(k,:,:) = cube(:,:,k)
    enddo
    deallocate(cube)
    
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
      


  function get_norm_psf_kernel(x,y,psf)
    implicit none
    real(fsp) :: get_norm_psf_kernel
    real(fsp), intent(in) :: x,y
    real(fsp), dimension(:,:) :: psf
    

    integer :: i,j,na,nb

    real(fsp) :: xnorm, ynorm
    real(fsp) :: sum

    na = size(psf,1)
    nb = size(psf,2)

    xnorm = (x-xoffset)/xscale
    ynorm = (y-yoffset)/yscale
    sum = 0._fsp
!$omp simd &
!$omp private(i,j) &
!$omp reduction(+:sum) &
!$omp collapse(2)
    do j=1,nb
       do i=1,na
          psf(i,j) = kernels(1,i,j) + xnorm * kernels(2,i,j) + ynorm * kernels(3,i,j) &
               + xnorm*ynorm * kernels(4,i,j) + xnorm*xnorm * kernels(5,i,j) + ynorm*ynorm * kernels(6,i,j)
          sum = sum + psf(i,j)
       enddo
    enddo
!$omp end simd          

    get_norm_psf_kernel = sum
    
  end function  get_norm_psf_kernel
  
  

end module psfeed
