!
! $Id: error_check.f90,v 1.9 2003/11/18 15:27:12 k-yukino Exp $
!
! SYNOPSIS
! 
!   call someting(arguments,name,error) 
!   call error_check(subroutine_name,error)
!
!    This subroutine checks your subroutines result. This is very useful to 
!   check it . I recommend for this to call after every subroutine was called .
!    First you set variable "name" for subroutine name that you want to 
!   check . After call subroutine that you want to check , you should call 
!   this error_check subroutine . 
!
! DESCRIPTION
!
!    If error<0 , then it displays name of subroutine that commit an error and
!   particular error number . And then this subroutine stops to execute this 
!   program . 
!    Else this subroutine sets error=0 and return to biginning of call .
!
subroutine error_check(name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  character(32),intent(in) :: name
  integer,intent(in) :: error
!
  if( error<0 ) then 
     write(*,*)
     write(*,'(a14,i5)') " (ToT)  MYRANK=",MYRANK
     write(*,'(a37,a32)') " (ToT)  an error occurred subroutine ",name
     write(*,'(a17,i5)') " (ToT)  error No.",error

     if( NUMBER_PE>=2 ) then
        call MPI_ABORT(MPI_COMM_WORLD,IERROR)
     end if

     stop
  end if

end subroutine error_check



#ifdef DEBUG_ERROR_CHECK
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,external :: iargc
  character(32) :: string           ! 引数を入れるバッファ
  character(32) :: name
  integer :: error
! init
  name="main" ; error=0
  IERROR=0

  call MPI_INIT(IERROR)

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERROR)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_PE,IERROR)

  name="TEST" ; error=-1

  call MPI_BCAST(name,32,MPI_CHARACTER,0,MPI_COMM_WORLD,IERROR) 
  call MPI_BCAST(error,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR) 

  call error_check(name,error)

  call MPI_FINALIZE(IERROR)

end program 
#endif /* DEBUG_ERROR_CHECK */
