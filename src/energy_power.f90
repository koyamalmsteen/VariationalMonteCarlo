!
! $Id: energy_power.f90,v 1.4 2004/02/15 10:52:35 k-yukino Exp ju448 $
!

!
! すでに集められた<alpha|を使って、パワーをかける。
!

!
! E_pow1=<Psi_1 |H|Psi_1 >/<Psi_1 |Psi_1 >
!  =C1^2 <phi_0 |H^3 |phi_0 >+2C1 <phi_0 |H^2 |phi_0 >+<phi_0 |H|phi_0 >
!   --------------------------------------------------------------------
!   C1^2 <phi_0 |H^2 |phi_0 >+2C1 <phi_0 |H|phi_0 >+<phi_0 |phi_0 >  
!

!
! E_pow2=<Psi_2 |H|Psi_2 >/<Psi_2 |Psi_2 >
!

#include "parameter.h"

subroutine energy_power(vector_up,vector_down,alpha_psi1,&
     alpha_h_psi1,h0,h1,h2,h3,h4,h5,alpha_h_h_psi1,&
     alpha_h_h_h_psi1,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1
  real(8),intent(out) :: h0,h1,h2,h3,h4,h5
  real(8),intent(out) :: alpha_h_h_psi1,alpha_h_h_h_psi1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  h0=0 ; h1=0 ; h2=0 ; h3=0 ; h4=0 ; h5=0
  name="energy_power" ; error=0
! init [local]
  alpha_h_h_psi1=0 ; alpha_h_h_h_psi1=0

!
! メインルーチン
!

  call calc_element_pow3(vector_up,vector_down,alpha_h_h_psi1,name,error)
  call error_check(name,error)

  if( POWER==2 ) then
! pow5
     call calc_element_pow5(vector_up,vector_down,&
          alpha_h_h_h_psi1,name,error)
     call error_check(name,error)
  end if
!
!
!
  h0=(alpha_psi1**2)/(alpha_psi1**2)
  h1=(alpha_psi1*alpha_h_psi1)/(alpha_psi1**2)
  h2=(alpha_h_psi1**2)/(alpha_psi1**2)
  h3=(alpha_h_psi1*alpha_h_h_psi1)/(alpha_psi1**2)

  if( POWER==2 ) then
     h4=(alpha_h_h_psi1**2)/(alpha_psi1**2)
     h5=(alpha_h_h_psi1*alpha_h_h_h_psi1)/(alpha_psi1**2)
  end if

end subroutine energy_power



#ifdef DEBUG_CALC_ELEMENT_POWER
program main
  use global_variables 
  implicit none
  
  include "mpif.h"
  
  real(8) :: result
  character(32) :: name
  integer :: error
  integer :: k_count,tmp,ier,IERROR
  integer,external :: count_number_xi
!
  name="main" ; error=0
  
  call MPI_INIT(IERROR)

  call init_distance()

  tmp=count_number_xi()

  allocate(vtable(tmp),stat=ier)
  call stat_check("vtable","main",1,ier)
  do k_count=1,tmp
     vtable(k_count)%way=0
     vtable(k_count)%xi=0
  end do

! set variation parameter
  do k_count=1,tmp

  end do

  call calc_element_power(result,name,error)  
  call error_check(name,error)

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ELEMENT_POWER */
