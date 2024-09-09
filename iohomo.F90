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

module iohomo
  use precision, only : fsp
  implicit none

  private
  
  integer, parameter :: lenarg = 80
  character(len=lenarg) :: inlist
  
  type fitsfile
     integer :: n
     character(len=lenarg), dimension(:), allocatable :: files
  end type fitsfile
  
  public fitsfile, get_arg_files, free_fitsfile
  
  
contains

  subroutine free_fitsfile(afits)
    implicit none
    type(fitsfile) :: afits

    if (allocated(afits%files)) deallocate(afits%files)
    
  end subroutine free_fitsfile
    

  subroutine get_arg_files(lists)
    implicit  none

    type(fitsfile), dimension(:), allocatable :: lists
    
    integer :: count
    character(len=lenarg) :: arg, buffer
    character(len=:), allocatable :: inlist
    logical :: isthere

    
    integer, parameter :: unit =110
    
    integer :: i,j
    integer :: ios
    
    count = command_argument_count()

    select case(count)

    case(3)

       allocate(lists(count))
       
       do i=1,count
          call get_command_argument(i,arg)

          if (index(arg,'@').eq.1) then
             if (allocated(inlist)) deallocate(inlist)
             inlist = arg(2:len_trim(arg))
             inquire(file=inlist,exist=isthere)
             if (.not.isthere) then
                write(*,*)'get_arg_files: files ',inlist,' not found!'
                stop
             endif
             
             open(unit,file=inlist,action='read',form='formatted')
             j=0
             ios=0
             do
                read(unit,iostat=ios,fmt='(A)') buffer
                if (ios.ne.0) exit
                j=j+1
             enddo
             rewind(unit)
             
             lists(i)%n=j
             allocate(lists(i)%files(j))
             do j=1,lists(i)%n
                read(unit,iostat=ios,fmt='(A)')lists(i)%files(j)
                if (ios.ne.0) stop 'get_arg_files: wrong counting!'
             enddo
             close(unit)
             
          else

             lists(i)%n=1
             allocate(lists(i)%files(1))
             lists(i)%files(1) = trim(arg)

          endif       
                                    
       enddo

       do i=1,count-1
          if (lists(i)%n .ne.lists(count)%n) then
             stop 'get_arg_files: input lists should be of same size!'
          endif
       enddo

       if (allocated(inlist)) deallocate(inlist)
       
    case default
       write(*,*)'Usage:    homogenix @inlistfiles @inkernelfiles @outlistfiles'
       write(*,*)'          homogenix inimage kernelcube outimage              '
       stop
    end select

        
  end subroutine get_arg_files



end module iohomo
