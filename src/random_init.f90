!
! $Id: random_init.f90,v 1.8 2003/11/18 08:05:41 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine random_init(ic,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(2),intent(in) :: ic
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="random_init" ; error=0
  
! デバッグ用
  call random_seed(put=(/ic/))

end subroutine random_init


#ifdef DEBUG_RANDOM_INIT
program main
  use global_variables
  implicit none

  include "mpif.h"

  real(8) :: ransuu
  character(32) :: name
  integer :: error
! MPI関連
  integer :: IERROR
! init
  ransuu=0 ; name="main" ; error=0
  
  call MPI_INIT(IERROR)
  
  call random_init()
  call fortran_random(ransuu,name,error) 
  
  write(*,*) "ransuu=",ransuu

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_RANDOM_INIT */
