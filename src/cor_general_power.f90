!
! $Id: cor_general_power.f90,v 1.3 2003/11/18 08:05:41 k-yukino Exp $
!

#include "parameter.h"

subroutine cor_general_power(alpha_psi0,alpha_psi1,alpha_h_psi1,&
     alpha_h_h_psi1,lambda,vector_up,vector_down,&
     tmp_sxsy_nume1,tmp_sxsy_nume2,tmp_sxsy_nume3_1,&
     tmp_sxsy_nume3_2,tmp_sxsy_nume4,tmp_sxsy_nume5,&
     tmp_super_nume1,tmp_super_nume2,tmp_super_nume3_1,&
     tmp_super_nume3_2,tmp_super_nume4,tmp_super_nume5,name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  real(8),intent(in) :: alpha_psi0,alpha_psi1,alpha_h_psi1,&
       alpha_h_h_psi1,lambda
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),dimension(max_number_xi),intent(in) :: tmp_sxsy_nume1,&
       tmp_sxsy_nume2,tmp_sxsy_nume3_1,&
       tmp_sxsy_nume3_2,tmp_sxsy_nume4,tmp_sxsy_nume5
  real(8),dimension(max_number_xi),intent(in) :: tmp_super_nume1,&
       tmp_super_nume2,tmp_super_nume3_1,&
       tmp_super_nume3_2,tmp_super_nume4,tmp_super_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="cor_power" ; error=0
! local init

!
! sxsy
!

  call sxsy_power(alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,lambda,&
       vector_up,vector_down,tmp_sxsy_nume1,tmp_sxsy_nume2,&
       tmp_sxsy_nume3_1,tmp_sxsy_nume3_2,tmp_sxsy_nume4,tmp_sxsy_nume5,&
       name,error)
  call error_check(name,error)

!
! super
!
  if( JIGEN==2 ) then
     call super_power(alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
          lambda,vector_up,vector_down,tmp_super_nume1,tmp_super_nume2,&
          tmp_super_nume3_1,tmp_super_nume3_2,tmp_super_nume4,&
          tmp_super_nume5,name,error)
     call error_check(name,error)
  end if

end subroutine cor_general_power
