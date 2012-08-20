!
! $Id: state_include.f90,v 1.8 2002/12/20 05:41:06 k-yukino Exp $
!

#include "parameter.h"

integer function state_include(total_sigma_electron,aaa,ket_vector_sigma)
  implicit none

  include "mpif.h"
  
  integer,intent(in) :: total_sigma_electron,aaa
  integer,dimension(total_sigma_electron),intent(in) :: ket_vector_sigma
! local
  integer :: i_count
! argument check
  if( total_sigma_electron<0 ) then
     write(*,*) "an error occured ,state_include ,No.1" ; stop
  else if( aaa<=0 ) then  
     write(*,*) "an error occured ,state_include ,No.2" ; stop
  end if

  do i_count=1,total_sigma_electron
     if( ket_vector_sigma(i_count)==aaa ) then
        state_include=1 ; return
     end if
  end do

! ket_vector_sigmaのいずれもaaaと等しくなければ
  state_include=0

end function state_include



#ifdef DEBUG_STATE_INCLUDE
program main
  implicit none

  include "mpif.h"

  integer :: total_sigma_electron,aaa 
  integer,external :: state_include
  integer,dimension(:),allocatable :: gamma_sigma
! local
  integer :: i_count,ier 
! MPI関連
  integer :: IERROR
! init
  i_count=0 ; ier=0
  IERROR=0

  call MPI_INIT(IERROR)

  write(*,*) "total_sigma_electron >" ; read(*,*) total_sigma_electron
  
  allocate(gamma_sigma(total_sigma_electron),stat=ier)
  call stat_check("gamma_sigma","state_include",1,ier)

  do i_count=1,total_sigma_electron
     write(*,*) "gamma_sigma >" ; read(*,*) gamma_sigma(i_count)
  end do

  write(*,*) "aaa >" ; read(*,*) aaa

  write(*,*) "result=",state_include(total_sigma_electron,aaa,gamma_sigma)
  
  call MPI_FINALIZE(IERROR)

end program main


#endif /* DEBUG_STATE_INCLUDE */
