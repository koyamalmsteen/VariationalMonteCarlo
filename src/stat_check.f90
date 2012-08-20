!
! $Id: stat_check.f90,v 1.11 2003/11/20 15:06:57 k-yukino Exp $
!

subroutine stat_check(variable_name,subroutine_name,status,ier)
  implicit none

  include "mpif.h"

  character(32),intent(in) :: variable_name
  character(32),intent(in) :: subroutine_name
  integer,intent(in) :: status              !       1:allocate , 2:deallocate
  integer,intent(in) :: ier
! MPI関連
  integer :: COMM,NUMBER_PE,MYRANK,IERROR

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERROR)

  COMM=MPI_COMM_WORLD
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_PE,IERROR)

! argument check
  if( ier<0 ) then 
     write(*,*) 
     write(*,'(a14,i5)') " (ToT) MYRANK=",MYRANK
     write(*,'(a56)') " (ToT) fatel error occured in stat_check(argument check)"
     
     if( NUMBER_PE>=2 ) then
        call MPI_ABORT(COMM,IERROR)
     end if

     stop
  else if( ier>0 ) then
     if( status==1 ) then
        write(*,*) 
        write(*,'(a14,i5)') " (ToT) MYRANK=",MYRANK
        write(*,'(a34,a32,a4,a32)') " (ToT) an error occurred allocate ",&
                             variable_name," in ",subroutine_name
        write(*,'(a16,i5)') " (ToT) error No.",ier
        write(*,*) 

        if( NUMBER_PE>=2 ) then
           call MPI_ABORT(COMM,IERROR)
        end if

        stop
     else if( status==2 ) then
        write(*,*) 
        write(*,'(a14,i5)') " (ToT) MYRANK=",MYRANK
        write(*,'(a37,a32,a4,a32)') " (ToT) an error occurred deallocate ",&
                                    variable_name," in ",subroutine_name
        write(*,'(a16,i5)') " (ToT) error No.",ier
        write(*,*)  
        
        if( NUMBER_PE>=2 ) then
           call MPI_ABORT(COMM,IERROR)
        end if

        stop
     else 
        write(*,*) 
        write(*,'(a14,i5)') " (ToT) MYRANK=",MYRANK
        write(*,'(a37)') " (ToT) an error occuered , Bad status" 
        
        if( NUMBER_PE>=2 ) then
           call MPI_ABORT(COMM,IERROR)
        end if

        stop
     end if
  else if( ier==0 ) then
     return
  end if

end subroutine stat_check



#ifdef DEBUG_STAT_CHECK
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer :: ier
  real(8),dimension(:),allocatable :: test 
! local  
! init
  IERROR=0

  call MPI_INIT(IERROR)

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERROR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMBER_PE,IERROR)
  
  allocate(test(10),stat=ier)
  call stat_check("test","main",1,ier)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_STAT_CHECK */
