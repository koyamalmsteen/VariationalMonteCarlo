!
! $Id: input_cor_general_power1.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine input_cor_general_power1(pow1_c0,pow1_c1,alpha_psi0,&
     alpha_psi1,alpha_h_psi1,&
     alpha_h_h_psi1,lambda,vector_up,vector_down,&
     sz_nume1,sz_nume2,sz_nume3_1,&
     sz_nume3_2,sz_nume4,sz_nume5,&
     sxsy_nume1,sxsy_nume2,sxsy_nume3_1,&
     sxsy_nume3_2,sxsy_nume4,sxsy_nume5,&
     super_nume1,super_nume2,super_nume3_1,&
     super_nume3_2,super_nume4,super_nume5,&
     name,error)
  use global_variables
  implicit none

  real(8),intent(in) :: pow1_c0,pow1_c1
  real(8),intent(in) :: alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
       lambda
  real(8),dimension(max_number_xi),intent(in) :: sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(max_number_xi),intent(in) :: sxsy_nume1,sxsy_nume2,sxsy_nume3_1,sxsy_nume3_2,sxsy_nume4,sxsy_nume5
  real(8),dimension(max_number_xi),intent(in) :: super_nume1,super_nume2,super_nume3_1,super_nume3_2,super_nume4,super_nume5
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down

  character(32),intent(out) :: name
  integer,intent(out) :: error
! init
  name="input_cor_general_power1" ; error=0

! sxsy
  result_sxsy1=result_sxsy1+( sxsy_nume3_1*(pow1_c1**2) &
       /(alpha_psi1**2) &
       +2*sxsy_nume2*(pow1_c0*pow1_c1) &
       /(alpha_psi1**2) &
       +sxsy_nume1*(pow1_c0**2) &
       /(alpha_psi1**2) ) 
! s
  result_s1=result_s1+( (sz_nume3_1+sxsy_nume3_1)*(pow1_c1**2) &
       /(alpha_psi1**2) &
       +2*(sz_nume2+sxsy_nume2)*(pow1_c0*pow1_c1) &
       /(alpha_psi1**2) &
       +(sz_nume1+sxsy_nume1)*(pow1_c0**2) &
       /(alpha_psi1**2) ) 
! super
  result_super1=result_super1+( super_nume3_1*(pow1_c1**2) &
       /(alpha_psi1**2) &
       +2*super_nume2*(pow1_c0*pow1_c1) &
       /(alpha_psi1**2) &
       +super_nume1*(pow1_c0**2) &
       /(alpha_psi1**2) ) 

end subroutine input_cor_general_power1
