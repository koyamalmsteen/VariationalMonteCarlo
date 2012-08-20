!
! $Id: input_splited_cor_general_power.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine input_splited_cor_general_power(sample_number,pow1_c0,pow1_c1,&
     pow2_c0,pow2_c1,pow2_c2,alpha_psi0,alpha_psi1,alpha_h_psi1,&
     alpha_h_h_psi1,lambda,vector_up,vector_down,&
     sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,&
     sxsy_nume1,sxsy_nume2,sxsy_nume3_1,sxsy_nume3_2,sxsy_nume4,sxsy_nume5,&
     super_nume1,super_nume2,super_nume3_1,super_nume3_2,&
     super_nume4,super_nume5,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  real(8),intent(in) :: pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2
  real(8),intent(in) :: alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
       lambda
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),dimension(max_number_xi),intent(in) :: sz_nume1,sz_nume2,sz_nume3_1,&
       sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(max_number_xi),intent(in) :: sxsy_nume1,sxsy_nume2,&
       sxsy_nume3_1,sxsy_nume3_2,sxsy_nume4,sxsy_nume5
  real(8),dimension(max_number_xi),intent(in) :: super_nume1,super_nume2,&
       super_nume3_1,super_nume3_2,super_nume4,super_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_splited_cor_general_power" ; error=0

  if( mod(sample_number,5)==1 ) then
     do i_count=1,max_number_xi
        splited_sxsy1(1,i_count)=splited_sxsy1(1,i_count) &
             +( sxsy_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sxsy_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sxsy_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 

        splited_s1(1,i_count)=splited_s1(1,i_count) &
             +( (sz_nume3_1(i_count)+sxsy_nume3_1(i_count))*(pow1_c1**2) &
             /(alpha_psi1**2) &
             +2*(sz_nume2(i_count)+sxsy_nume2(i_count))*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +(sz_nume1(i_count)+sxsy_nume1(i_count))*(pow1_c0**2) &
             /(alpha_psi1**2) ) 

        splited_super1(1,i_count)=splited_super1(1,i_count) &
             +( super_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*super_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +super_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 

        if( POWER==2 ) then
           splited_sxsy2(1,i_count)=splited_sxsy2(1,i_count) &
                +( sxsy_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sxsy_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sxsy_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sxsy_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sxsy_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sxsy_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_s2(1,i_count)=splited_s2(1,i_count) &
                +(sz_nume5(i_count)+sxsy_nume5(i_count))*(pow2_c2**2) &
                /(alpha_psi1**2) &
                +2*(sz_nume4(i_count)+sxsy_nume4(i_count)) &
                *(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*(sz_nume3_2(i_count)+sxsy_nume3_2(i_count)) &
                *(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +(sz_nume3_1(i_count)+sxsy_nume3_1(i_count)) &
                *(pow2_c1**2)/(alpha_psi1**2) &
                +2*(sz_nume2(i_count)+sxsy_nume2(i_count)) &
                *(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +(sz_nume1(i_count)+sxsy_nume1(i_count)) &
                *(pow2_c0**2)/(alpha_psi1**2)
           splited_super2(1,i_count)=splited_super2(1,i_count) &
                +( super_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*super_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*super_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +super_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*super_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +super_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end if
     end do
  else if( mod(sample_number,5)==2 ) then
     do i_count=1,max_number_xi
        splited_sxsy1(2,i_count)=splited_sxsy1(2,i_count) &
             +( sxsy_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sxsy_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sxsy_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        splited_s1(2,i_count)=splited_s1(2,i_count) &
             +( (sz_nume3_1(i_count)+sxsy_nume3_1(i_count))*(pow1_c1**2) &
             /(alpha_psi1**2) &
             +2*(sz_nume2(i_count)+sxsy_nume2(i_count))*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +(sz_nume1(i_count)+sxsy_nume1(i_count))*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        splited_super1(2,i_count)=splited_super1(2,i_count) &
             +( super_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*super_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +super_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        if( POWER==2 ) then
           splited_sxsy2(2,i_count)=splited_sxsy2(2,i_count) &
                +( sxsy_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sxsy_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sxsy_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sxsy_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sxsy_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sxsy_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_s2(2,i_count)=splited_s2(2,i_count) &
                +(sz_nume5(i_count)+sxsy_nume5(i_count)) &
                *(pow2_c2**2)/(alpha_psi1**2) &
                +2*(sz_nume4(i_count)+sxsy_nume4(i_count)) &
                *(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*(sz_nume3_2(i_count)+sxsy_nume3_2(i_count)) &
                *(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +(sz_nume3_1(i_count)+sxsy_nume3_1(i_count)) &
                *(pow2_c1**2)/(alpha_psi1**2) &
                +2*(sz_nume2(i_count)+sxsy_nume2(i_count)) &
                *(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +(sz_nume1(i_count)+sxsy_nume1(i_count)) &
                *(pow2_c0**2)/(alpha_psi1**2)
           splited_super2(2,i_count)=splited_super2(2,i_count) &
                +( super_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*super_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*super_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +super_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*super_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +super_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end if
     end do
  else if( mod(sample_number,5)==3 ) then
     do i_count=1,max_number_xi
        splited_sxsy1(3,i_count)=splited_sxsy1(3,i_count) &
             +( sxsy_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sxsy_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sxsy_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        splited_s1(3,i_count)=splited_s1(3,i_count) &
             +( (sz_nume3_1(i_count)+sxsy_nume3_1(i_count))*(pow1_c1**2) &
             /(alpha_psi1**2) &
             +2*(sz_nume2(i_count)+sxsy_nume2(i_count))*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +(sz_nume1(i_count)+sxsy_nume1(i_count))*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        splited_super1(3,i_count)=splited_super1(3,i_count) &
             +( super_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*super_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +super_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        if( POWER==2 ) then
           splited_sxsy2(3,i_count)=splited_sxsy2(3,i_count) &
                +( sxsy_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sxsy_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sxsy_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sxsy_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sxsy_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sxsy_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_s2(3,i_count)=splited_s2(3,i_count) &
                +(sz_nume5(i_count)+sxsy_nume5(i_count)) &
                *(pow2_c2**2)/(alpha_psi1**2) &
                +2*(sz_nume4(i_count)+sxsy_nume4(i_count)) &
                *(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*(sz_nume3_2(i_count)+sxsy_nume3_2(i_count)) &
                *(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +(sz_nume3_1(i_count)+sxsy_nume3_1(i_count)) &
                *(pow2_c1**2)/(alpha_psi1**2) &
                +2*(sz_nume2(i_count)+sxsy_nume2(i_count)) &
                *(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +(sz_nume1(i_count)+sxsy_nume1(i_count)) &
                *(pow2_c0**2)/(alpha_psi1**2)
           splited_super2(4,i_count)=splited_super2(4,i_count) &
                +( super_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*super_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*super_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +super_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*super_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +super_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end if
     end do
  else if( mod(sample_number,5)==4 ) then
     do i_count=1,max_number_xi
        splited_sxsy1(5,i_count)=splited_sxsy1(5,i_count) &
             +( sxsy_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sxsy_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sxsy_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        splited_s1(5,i_count)=splited_s1(5,i_count) &
             +( (sz_nume3_1(i_count)+sxsy_nume3_1(i_count))*(pow1_c1**2) &
             /(alpha_psi1**2) &
             +2*(sz_nume2(i_count)+sxsy_nume2(i_count))*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +(sz_nume1(i_count)+sxsy_nume1(i_count))*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        splited_super1(5,i_count)=splited_super1(5,i_count) &
             +( super_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*super_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +super_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) ) 
        if( POWER==2 ) then
           splited_sxsy2(5,i_count)=splited_sxsy2(5,i_count) &
                +( sxsy_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sxsy_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sxsy_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sxsy_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sxsy_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sxsy_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_s2(5,i_count)=splited_s2(5,i_count) &
                +(sz_nume5(i_count)+sxsy_nume5(i_count)) &
                *(pow2_c2**2)/(alpha_psi1**2) &
                +2*(sz_nume4(i_count)+sxsy_nume4(i_count)) &
                *(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*(sz_nume3_2(i_count)+sxsy_nume3_2(i_count)) &
                *(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +(sz_nume3_1(i_count)+sxsy_nume3_1(i_count)) &
                *(pow2_c1**2)/(alpha_psi1**2) &
                +2*(sz_nume2(i_count)+sxsy_nume2(i_count)) &
                *(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +(sz_nume1(i_count)+sxsy_nume1(i_count)) &
                *(pow2_c0**2)/(alpha_psi1**2)
           splited_super2(5,i_count)=splited_super2(5,i_count) &
                +( super_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*super_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*super_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +super_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*super_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +super_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end if
     end do
  else if( mod(sample_number,5)==0 ) then
     do i_count=1,max_number_xi
        result_sxsy1=result_sxsy1+( sxsy_nume3_1*(pow1_c1**2) &
             /(alpha_psi1**2)+2*sxsy_nume2*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2)+sxsy_nume1*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        result_s1=result_s1+( (sz_nume3_1+sxsy_nume3_1)*(pow1_c1**2) &
             /(alpha_psi1**2)+2*(sz_nume2+sxsy_nume2)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2)+(sz_nume1+sxsy_nume1)*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        result_super1=result_super1+( super_nume3_1*(pow1_c1**2) &
             /(alpha_psi1**2)+2*super_nume2*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2)+super_nume1*(pow1_c0**2) &
             /(alpha_psi1**2) ) 
        if( POWER==2 ) then
           splited_sxsy2(1,i_count)=splited_sxsy2(1,i_count) &
                +( sxsy_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sxsy_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sxsy_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sxsy_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sxsy_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sxsy_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_s2(1,i_count)=splited_s2(1,i_count) &
                +(sz_nume5(i_count)+sxsy_nume5(i_count)) &
                *(pow2_c2**2)/(alpha_psi1**2) &
                +2*(sz_nume4(i_count)+sxsy_nume4(i_count)) &
                *(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*(sz_nume3_2(i_count)+sxsy_nume3_2(i_count)) &
                *(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +(sz_nume3_1(i_count)+sxsy_nume3_1(i_count)) &
                *(pow2_c1**2)/(alpha_psi1**2) &
                +2*(sz_nume2(i_count)+sxsy_nume2(i_count)) &
                *(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +(sz_nume1(i_count)+sxsy_nume1(i_count)) &
                *(pow2_c0**2)/(alpha_psi1**2)
           splited_super2(1,i_count)=splited_super2(1,i_count) &
                +( super_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*super_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*super_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +super_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*super_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +super_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end if
     end do
  end if

end subroutine input_splited_cor_general_power
