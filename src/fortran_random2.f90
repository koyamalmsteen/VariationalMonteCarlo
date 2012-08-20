!
! $Id: fortran_random2.f90,v 1.3 2003/02/28 01:59:29 k-yukino Exp $
!
subroutine fortran_random2(aaa,bbb,ransuu_int,name,error)
  implicit none

  include "mpif.h"

  integer,intent(in) :: aaa,bbb
  integer,intent(out) :: ransuu_int
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: ransuu
! init
  ransuu_int=0 ; name="fortran_random2" ; error=0
  ransuu=0

  if( aaa<bbb ) then
     call fortran_random(ransuu,name,error)
     call error_check(name,error)
     ransuu_int=aaa+int((bbb-aaa+1)*ransuu)
  else if( aaa>bbb ) then
     call fortran_random(ransuu,name,error)
     call error_check(name,error)
     ransuu_int=bbb+int((aaa-bbb+1)*ransuu)
  else
     ransuu_int=aaa
  end if

end subroutine fortran_random2



#ifdef DEBUG_FORTRAN_RANDOM2
program main
  implicit none

  include "mpif.h"
  
  character(32) :: name
  integer :: error
! local
  integer :: ransuu_int
! MPI関連
  integer :: IERROR,i_count
! init  
  name="main" ; error=0
  ransuu_int=0

  call MPI_INIT(IERROR)

  call random_init(name,error)
  call error_check(name,error)

  call fortran_random2(1,10,ransuu_int,name,error)
  call error_check(name,error)

  write(*,*) "ransuu_int=",ransuu_int

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_FORTRAN_RANDOM2 */
