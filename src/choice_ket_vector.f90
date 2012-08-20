!
! $Id: choice_ket_vector.f90,v 1.13 2003/11/05 13:32:46 k-yukino Exp $
!
#include "parameter.h"

subroutine choice_ket_vector(gamma_up,gamma_down,number_ket_vector_up,&
                             number_ket_vector_down,available_ket_up_table,&
                             available_ket_down_table,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down
  integer,intent(out) :: number_ket_vector_up,number_ket_vector_down
  integer,dimension(2*JIGEN*TOTAL_UP_ELECTRON+1,TOTAL_UP_ELECTRON),intent(out) :: available_ket_up_table
  integer,dimension(2*JIGEN*TOTAL_DOWN_ELECTRON+1,TOTAL_DOWN_ELECTRON),intent(out) :: available_ket_down_table
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier,i_count,j_count,k_count,l_count,m_count
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
! init
  number_ket_vector_up=0 ; number_ket_vector_down=0
  available_ket_up_table=0 ; available_ket_down_table=0
  name="choice_ket_vector" ; error=0
! init [local]
  ier=0 ; i_count=0 ; j_count=0 ; k_count=0 ; l_count=0 ; m_count=0
  site_table=0

! $B",%9%T%s(B($B%5%$%HI=<($K$9$k(B)

  site_table=0
  do i_count=1,TOTAL_UP_ELECTRON
     site_table(gamma_up(i_count))=1
  end do

  l_count=1                              ! $B%V%i%Y%/%H%k<+?H$rF~$l$k$N$G(B 

  do i_count=1,TOTAL_UP_ELECTRON ! $B0\F085(B
     do j_count=1,TOTAL_SITE_NUMBER ! $B0\F0@h(B
        if( neighbor_table(gamma_up(i_count),j_count)==1 ) then! $B"+NY@\$J$i0\F02D(B
           if( site_table(j_count)==0 ) then ! $B$b$70\F0@h$,6u$J$i(B
              do k_count=1,TOTAL_UP_ELECTRON
                 available_ket_up_table(l_count,k_count)=gamma_up(k_count)
              end do
              available_ket_up_table(l_count,i_count)=j_count
              l_count=l_count+1
           end if
        end if
     end do
  end do

  do i_count=1,TOTAL_UP_ELECTRON  ! $B%V%i%Y%/%H%k<+?H(B
     available_ket_up_table(l_count,i_count)=gamma_up(i_count)
  end do

  number_ket_vector_up=l_count

! $B"-%9%T%s(B
  site_table=0
  do i_count=1,TOTAL_DOWN_ELECTRON
     site_table(gamma_down(i_count))=1
  end do

  l_count=1                              ! $B%V%i%Y%/%H%k<+?H$rF~$l$k$N$G(B 

  do i_count=1,TOTAL_DOWN_ELECTRON ! $B0\F085(B
     do j_count=1,TOTAL_SITE_NUMBER ! $B0\F0@h(B
        if( neighbor_table(gamma_down(i_count),j_count)==1 ) then! $B"+NY@\$J$i0\F02D(B
           if( site_table(j_count)==0 ) then ! $B$b$70\F0@h$,6u$J$i(B
              do k_count=1,TOTAL_DOWN_ELECTRON
                 available_ket_down_table(l_count,k_count)=gamma_down(k_count)
              end do
              available_ket_down_table(l_count,i_count)=j_count
              l_count=l_count+1
           end if
        end if
     end do
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON  ! $B%V%i%Y%/%H%k<+?H(B
     available_ket_down_table(l_count,i_count)=gamma_down(i_count)
  end do

  number_ket_vector_down=l_count

end subroutine choice_ket_vector



#ifdef DEBUG_CHOICE_KET_VECTOR
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer :: number_ket_vector_up,number_ket_vector_down
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer,dimension(2*JIGEN*TOTAL_UP_ELECTRON+1,TOTAL_UP_ELECTRON) :: available_ket_up_table
  integer,dimension(2*JIGEN*TOTAL_DOWN_ELECTRON+1,TOTAL_DOWN_ELECTRON) :: available_ket_down_table
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count,j_count
! MPI$B4XO"(B
  integer :: IERROR
! init
  name="main" ; error=0 ; ier=0 ; i_count=0 ; j_count=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "gamma_up >" ; read(*,*) gamma_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "gamma_down >" ; read(*,*) gamma_down(i_count)
  end do

  call choice_ket_vector(gamma_up,gamma_down,number_ket_vector_up,&
                         number_ket_vector_down,available_ket_up_table,&
                         available_ket_down_table,name,error)
  call error_check(name,error)

  write(*,*) "available_ket_up_table=",available_ket_up_table
  write(*,*) "available_ket_down_table=",available_ket_down_table

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHOICE_KET_VECTOR */
