!
! $Id: fortran_random.f90,v 1.12 2003/06/04 10:07:47 k-yukino Exp $
!
#include "parameter.h"

subroutine fortran_random(ransuu,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(out) :: ransuu
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  ransuu=0 ; name="fortran_random" ; error=0

! main
  if( RANDOM_ALGORITHM==0 ) then
     call random_number(ransuu)
  else
     error=-11 ; return
  end if

end subroutine fortran_random



#ifdef DEBUG_FORTRAN_RANDOM
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

  call fortran_random(ransuu,name,error)
  call error_check(name,error)

  write(*,*) ransuu 

  call MPI_FINALIZE(IERROR)

end program main

#endif /* DEBUG_FORTRAN_RANDOM */
