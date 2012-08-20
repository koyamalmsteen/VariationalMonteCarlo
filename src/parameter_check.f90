!
! $Id: parameter_check.f90,v 1.19 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine parameter_check(name,error)
  use global_variables
  implicit none

  include "mpif.h"

  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="parameter_check" ; error=0

  if( JIGEN/=1 .and. JIGEN/=2 ) then
     error=-91 ; return
  else if( TOTAL_SITE_NUMBER<=0 ) then
     error=-92 ; return
  else if( TOTAL_UP_ELECTRON<0 ) then
     error=-93 ; return
  else if( TOTAL_SITE_NUMBER<TOTAL_UP_ELECTRON ) then
     error=-923 ; return
  else if( TOTAL_DOWN_ELECTRON<0 ) then
    error=-94 ; return
  else if( TOTAL_SITE_NUMBER<TOTAL_DOWN_ELECTRON ) then
     error=-924 ; return
  else if( TRANSFER>=0 ) then              ! 0も計算不可能
     error=-95 ; return
  else if( COULOMB<0 ) then
     error=-96 ; return
  else if( INIT_WAVE_VECTOR/=0 .and. INIT_WAVE_VECTOR/=1 .and. INIT_WAVE_VECTOR/=2 ) then
     error=-97 ; return
  else if( GAMMA_ALGORITHM/=0 .and. GAMMA_ALGORITHM/=1 ) then
     error=-98 ; return
  else if( PERIODIC/=0 .and. PERIODIC/=1 .and. PERIODIC/=2 ) then
     error=-99 ; return
!
!
  else if( RANDOM_ALGORITHM/=0 ) then            ! 今のところはこれだけ
     
  else if( HF_INITIAL_RHO/=0 .and. HF_INITIAL_RHO/=1 .and. &
          &HF_INITIAL_RHO/=2 .and. HF_INITIAL_RHO/=3 ) then
     error=-910 ; return
  else if( TRASH<0 ) then
     error=-911 ; return
  else if( MAX_HF_ITERATION<=0 ) then
     error=-912 ; return
  else if( HF_CONVERGENCE_CONDITION<=0 ) then
     error=-913 ; return
!
!
!
  else if( COR/=0 .and. COR/=1 .and. COR/=2 ) then
     error=-914 ; return
  else if( POWER/=0 .and. POWER/=1 .and. POWER/=2 ) then
     error=-915 ; return
  else if( NUMBER_XI<0 ) then
     error=-918 ; return
  else if( MAX_MONTECARLO_SAMPLE<=0 ) then
     error=-919 ; return
!
!
!
  else if( ELECTRIC_FIELD<0 ) then
     error=-921 ; return
  else if( FLAT_BAND/=0 .and. FLAT_BAND/=1 .and. FLAT_BAND/=2 &
       .and. FLAT_BAND/=3) then
     error=-925 ; return
  else if( (FLAT_BAND==1 .and. JIGEN/=2) .and. &
       (FLAT_BAND==2 .and. JIGEN/=2) .and.  &
       (FLAT_BAND==3 .and. JIGEN/=2) ) then
     error=-9251 ; return
  else if( (FLAT_BAND==1 .and. PERIODIC/=1) .and. &
       (FLAT_BAND==2 .and. PERIODIC/=1) .and. &
       (FLAT_BAND==2 .and. PERIODIC/=1) ) then
     error=-9252 ; return
  end if

end subroutine parameter_check



#ifdef DEBUG_PARAMETER_CHECK
program main
  implicit none

  include "mpif.h"

  character(32) :: name
  integer :: error
! MPI関連
  integer :: IERROR
! init
  name="main" ; error=0

  call MPI_INIT(IERROR)

! main

  call parameter_check(name,error)
  call error_check(name,error)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_PARAMETER_CHECK */
