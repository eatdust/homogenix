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

module convolutor
  use precision, only : fsp
  use psfeed, only : get_norm_psf_kernel, get_psf_xsize, get_psf_ysize
  implicit none

  private

  abstract interface
     function get_kernel(x,y,ker)
       use precision, only : fsp
       implicit none
       real(fsp) :: get_kernel
       real(fsp), intent(in) :: x,y
       real(fsp), dimension(:,:) :: ker
     end function  get_kernel

     function get_xsize()
       implicit none
       integer :: get_xsize
     end function get_xsize

     function get_ysize()
       implicit none
       integer :: get_ysize
     end function get_ysize
     
  end interface

  procedure(get_kernel), pointer :: ptr_get_kernel => get_norm_psf_kernel
  procedure(get_xsize), pointer :: ptr_get_xsize => get_psf_xsize
  procedure(get_ysize), pointer :: ptr_get_ysize => get_psf_ysize
  
  public pix_convolve

  

contains

  subroutine pix_convolve(image,omage)
    implicit none
    real(fsp), dimension(:,:), intent(in) :: image
    real(fsp), dimension(:,:), intent(out) :: omage

    integer :: na,nb,np,nq,nx,ny
    integer :: i,j,p,q
    integer :: imp,jmq

    real(fsp) :: x,y,norm
    real(fsp), dimension(:,:), allocatable :: psf
    real(fsp) :: outpix
   
    
    nx = size(image,1)
    ny = size(image,2)
    
    na = ptr_get_xsize()
    nb = ptr_get_ysize()
    
    np = (na-1)/2
    nq = (nb-1)/2
    
!$omp parallel &    
!$omp default(shared) &
!$omp private(j,y,i,x,q,p,imp,jmq) &
!$omp private(psf,norm,outpix)
    allocate(psf(na,nb))
!$omp do
    do j=1,ny
     y = real(j,fsp)
     do i=1,nx
        x = real(i,fsp)
        norm = ptr_get_kernel(x,y,psf)
        outpix = 0._fsp
        do q=-nq,nq
           jmq = max(1,min(j-q,ny))
!$omp simd &
!$omp private(p,imp) &
!$omp reduction(+:outpix)
           do p=-np,np
              imp = max(1,min(i-p,nx))
#ifndef NONORM
              outpix = psf(p+np+1,q+nq+1)*image(imp,jmq)/norm + outpix
#else
              outpix = psf(p+np+1,q+nq+1)*image(imp,jmq) + outpix
#endif              
           enddo
!$omp end simd
        enddo
        omage(i,j) = outpix
     enddo
  enddo
!$omp end do
  deallocate(psf)
!$omp end parallel
  
end subroutine pix_convolve


end module convolutor
