!
! $Id: make_fock_matrix.f90,v 1.15 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine make_fock_matrix(rho_up,rho_down,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),dimension(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),intent(in) :: &
                                                               rho_up,rho_down
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
! init
  unitary_d_up=0 ; unitary_d_down=0
  name="make_fock_matrix" ; error=0 ; i_count=0 ; j_count=0 

! argument check

! main
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=dble(TRANSFER)*neighbor_table(i_count,j_count)
        unitary_d_down(i_count,j_count)=dble(TRANSFER)*neighbor_table(i_count,j_count)
     end do
  end do

  do i_count=1,TOTAL_SITE_NUMBER
     unitary_d_up(i_count,i_count)=unitary_d_up(i_count,i_count) &
                                   +dble(COULOMB)*rho_down(i_count,i_count)
     unitary_d_down(i_count,i_count)=unitary_d_down(i_count,i_count) &
                                   +dble(COULOMB)*rho_up(i_count,i_count)
  end do

end subroutine make_fock_matrix



#ifdef DEBUG_MAKE_FOCK_MATRIX
program main
  use global_variables
  implicit none

  real(8),dimension(TOTAL_SITE_NUMBER) :: rho_up,rho_down
  character(32) :: name
  integer :: error
! local
  real(8) :: ransuu
  integer :: i_count,j_count
! MPI関連
  integer :: IERROR
! init
  unitary_d_up=0 ; unitary_d_down=0
  name="main" ; error=0 ; i_count=0 ; j_count=0 ; ransuu=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

! まずは適当に乱数を発生させ、初期のrhoを決定する。
  call random_init(name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_SITE_NUMBER
     call fortran_random(ransuu,name,error)
     call error_check(name,error) 
     rho_up(i_count)=ransuu

     call fortran_random(ransuu,name,error)
     call error_check(name,error) 
     rho_down(i_count)=ransuu
  end do

  call make_fock_matrix(rho_up,rho_down,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        write(*,*) "(fock matrix)unitary_d_up(",i_count,j_count,")=",&
                   unitary_d_up(i_count,j_count)
     end do
  end do

  write(*,*) "rho_down=",rho_down
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON
  write(*,*) "TRANSFER=",TRANSFER
  write(*,*) "COULOMB=",COULOMB

  call MPI_FINALIZE(IERROR)

end program main
#endif
