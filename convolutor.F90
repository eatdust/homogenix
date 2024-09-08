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

module convolutor
  use precision, only : fsp
  use psfeed, only : get_psf_kernel, get_psf_xsize, get_psf_ysize
  implicit none

  private

  abstract interface
     subroutine get_kernel(x,y,ker)
       use precision, only : fsp
       implicit none
       real(fsp), intent(in) :: x,y
       real(fsp), dimension(:,:) :: ker
     end subroutine get_kernel

     function get_xsize()
       implicit none
       integer :: get_xsize
     end function get_xsize

     function get_ysize()
       implicit none
       integer :: get_ysize
     end function get_ysize
     
  end interface

  procedure(get_kernel), pointer :: ptr_get_kernel => get_psf_kernel
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
    
    nx = size(image,1)
    ny = size(image,2)
    
    na = ptr_get_xsize()
    nb = ptr_get_ysize()
    
    np = (na-1)/2
    nq = (nb-1)/2
    
!$omp parallel &    
!$omp default(shared) &
!$omp private(j,y,i,x,q,p,imp,jmq) &
!$omp private(psf,norm)
    allocate(psf(na,nb))
!$omp do
    do j=1,ny
     y = real(j,fsp)
     do i=1,nx
        x = real(i,fsp)
        call ptr_get_kernel(x,y,psf)
        norm = sum(psf)
        omage(i,j) = 0._fsp
        do q=-nq,nq
           jmq = max(1,min(j-q,ny))
           do p=-np,np
              imp = max(1,min(i-p,nx))
              omage(i,j) = psf(p+np+1,q+nq+1)*image(imp,jmq)/norm + omage(i,j)
           enddo
        enddo
     enddo
  enddo
!$omp end do
  deallocate(psf)
!$omp end parallel
  
end subroutine pix_convolve


end module convolutor
