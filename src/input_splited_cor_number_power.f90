!
! $Id: input_splited_cor_number_power.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine input_splited_cor_number_power(sample_number,alpha_psi1,alpha_h_psi1,&
     alpha_h_h_psi1,pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2,&
     sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5,&
     charge_nume1,charge_nume2,charge_nume3_1,charge_nume3_2,charge_nume4,&
     charge_nume5,charge_density_nume1,charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2,charge_density_nume4,charge_density_nume5,&
     name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  real(8),intent(in) :: alpha_psi1,alpha_h_psi1,alpha_h_h_psi1
  real(8),intent(in) :: pow1_c0,pow1_c1,pow2_c0,pow2_c1,pow2_c2 
  real(8),dimension(max_number_xi),intent(in) :: sz_nume1,sz_nume2,sz_nume3_1,sz_nume3_2,sz_nume4,sz_nume5
  real(8),dimension(max_number_xi),intent(in) :: charge_nume1,charge_nume2,charge_nume3_1,charge_nume3_2,charge_nume4,charge_nume5
  real(8),dimension(TOTAL_SITE_NUMBER),intent(in) :: charge_density_nume1,&
       charge_density_nume2,charge_density_nume3_1,charge_density_nume3_2,&
       charge_density_nume4,charge_density_nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_splited_cor_number" ; error=0

  if( mod(sample_number,5)==1 ) then
     do i_count=1,max_number_xi
        splited_sz1(1,i_count)=splited_sz1(1,i_count) &
             +( sz_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sz_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sz_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )

        splited_charge1(1,i_count)=splited_charge1(1,i_count) &
             +( charge_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +charge_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
     end do
     do i_count=1,TOTAL_SITE_NUMBER

        splited_charge_density1(1,i_count)=splited_charge_density1(1,i_count) &
             +( charge_density_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_density_nume2(i_count)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +charge_density_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        splited_pole1(1)=splited_pole1(1) &
             +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
             *(pow1_c1**2) /(alpha_psi1**2) &
             +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
             *(pow1_c0*pow1_c1) /(alpha_psi1**2) &
             +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
             *(pow1_c0**2) /(alpha_psi1**2)
     end do

     if( POWER==2 ) then
        do i_count=1,max_number_xi
           splited_sz2(1,i_count)=splited_sz2(1,i_count) &
                +( sz_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sz_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sz_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sz_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sz_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sz_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_charge2(1,i_count)=splited_charge2(1,i_count) &
                +( charge_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
                +charge_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2 ) &
                +2*charge_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
                +charge_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end do

        do i_count=1,TOTAL_SITE_NUMBER
           splited_charge_density2(1,i_count)= &
                splited_charge_density2(1,i_count) &
                +( charge_density_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_density_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_density_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
                +charge_density_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*charge_density_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +charge_density_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_pole2(1)=splited_pole2(1) &
                +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count))&
                *(pow2_c2**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count))&
                *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count))&
                *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
                +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
                *(pow2_c1**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
                *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
                +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count))&
                *(pow2_c0**2) /(alpha_psi1**2) 
        end do
     end if

  else if( mod(sample_number,5)==2 ) then
     do i_count=1,max_number_xi

        splited_sz1(2,i_count)=splited_sz1(2,i_count) &
             +( sz_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sz_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sz_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )

        splited_charge1(2,i_count)=splited_charge1(2,i_count) &
             +( charge_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +charge_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
     end do

     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(2,i_count)=splited_charge_density1(2,i_count) &
             +( charge_density_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_density_nume2(i_count)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +charge_density_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        splited_pole1(2)=splited_pole1(2) &
             +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
             *(pow1_c1**2) /(alpha_psi1**2) &
             +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
             *(pow1_c0*pow1_c1) /(alpha_psi1**2) &
             +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
             *(pow1_c0**2) /(alpha_psi1**2)
     end do

     if( POWER==2 ) then
        do i_count=1,max_number_xi
           splited_sz2(2,i_count)=splited_sz2(2,i_count) &
                +( sz_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sz_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sz_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sz_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sz_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sz_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_charge2(2,i_count)=splited_charge2(2,i_count) &
                +( charge_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
                +charge_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2 ) &
                +2*charge_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
                +charge_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end do

        do i_count=1,TOTAL_SITE_NUMBER
           splited_charge_density2(2,i_count)= &
                splited_charge_density2(2,i_count) &
                +( charge_density_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_density_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_density_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
                +charge_density_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*charge_density_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +charge_density_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_pole2(2)=splited_pole2(2) &
                +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count))&
                *(pow2_c2**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count))&
                *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count))&
                *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
                +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
                *(pow2_c1**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
                *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
                +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count))&
                *(pow2_c0**2) /(alpha_psi1**2) 
        end do
     end if

  else if( mod(sample_number,5)==3 ) then
     do i_count=1,max_number_xi

        splited_sz1(3,i_count)=splited_sz1(3,i_count) &
             +( sz_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sz_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sz_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )

        splited_charge1(3,i_count)=splited_charge1(3,i_count) &
             +( charge_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +charge_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
     end do

     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(3,i_count)=splited_charge_density1(3,i_count) &
             +( charge_density_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_density_nume2(i_count)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +charge_density_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        splited_pole1(3)=splited_pole1(3) &
             +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
             *(pow1_c1**2) /(alpha_psi1**2) &
             +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
             *(pow1_c0*pow1_c1) /(alpha_psi1**2) &
             +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
             *(pow1_c0**2) /(alpha_psi1**2)
     end do

     if( POWER==2 ) then
        do i_count=1,max_number_xi
           splited_sz2(3,i_count)=splited_sz2(3,i_count) &
                +( sz_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sz_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sz_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sz_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sz_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sz_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )

           splited_charge2(3,i_count)=splited_charge2(3,i_count) &
                +( charge_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
                +charge_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2 ) &
                +2*charge_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
                +charge_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end do

        do i_count=1,TOTAL_SITE_NUMBER
           splited_charge_density2(3,i_count)= &
                splited_charge_density2(3,i_count) &
                +( charge_density_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_density_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_density_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
                +charge_density_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*charge_density_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +charge_density_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_pole2(3)=splited_pole2(3) &
                +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count))&
                *(pow2_c2**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count))&
                *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count))&
                *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
                +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
                *(pow2_c1**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
                *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
                +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count))&
                *(pow2_c0**2) /(alpha_psi1**2) 
        end do
     end if

  else if( mod(sample_number,5)==4 ) then
     do i_count=1,max_number_xi

        splited_sz1(4,i_count)=splited_sz1(4,i_count) &
             +( sz_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sz_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sz_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )

        splited_charge1(4,i_count)=splited_charge1(4,i_count) &
             +( charge_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +charge_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
     end do

     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(4,i_count)=splited_charge_density1(4,i_count) &
             +( charge_density_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_density_nume2(i_count)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +charge_density_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        splited_pole1(4)=splited_pole1(4) &
             +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
             *(pow1_c1**2) /(alpha_psi1**2) &
             +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
             *(pow1_c0*pow1_c1) /(alpha_psi1**2) &
             +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
             *(pow1_c0**2) /(alpha_psi1**2)
     end do

     if( POWER==2 ) then
        do i_count=1,max_number_xi
           splited_sz2(4,i_count)=splited_sz2(4,i_count) &
                +( sz_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sz_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sz_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sz_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sz_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sz_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )

           splited_charge2(4,i_count)=splited_charge2(4,i_count) &
                +( charge_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
                +charge_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2 ) &
                +2*charge_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
                +charge_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end do

        do i_count=1,TOTAL_SITE_NUMBER
           splited_charge_density2(4,i_count)= &
                splited_charge_density2(4,i_count) &
                +( charge_density_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_density_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_density_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
                +charge_density_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*charge_density_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +charge_density_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_pole2(4)=splited_pole2(4) &
                +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count))&
                *(pow2_c2**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count))&
                *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count))&
                *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
                +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
                *(pow2_c1**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
                *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
                +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count))&
                *(pow2_c0**2) /(alpha_psi1**2) 
        end do
     end if

  else if( mod(sample_number,5)==0 ) then
     do i_count=1,max_number_xi

        splited_sz1(5,i_count)=splited_sz1(5,i_count) &
             +( sz_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*sz_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +sz_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        
        splited_charge1(5,i_count)=splited_charge1(5,i_count) &
             +( charge_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_nume2(i_count)*(pow1_c0*pow1_c1)/(alpha_psi1**2) &
             +charge_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
     end do

     do i_count=1,TOTAL_SITE_NUMBER
        splited_charge_density1(5,i_count)=splited_charge_density1(5,i_count) &
             +( charge_density_nume3_1(i_count)*(pow1_c1**2)/(alpha_psi1**2) &
             +2*charge_density_nume2(i_count)*(pow1_c0*pow1_c1) &
             /(alpha_psi1**2) &
             +charge_density_nume1(i_count)*(pow1_c0**2)/(alpha_psi1**2) )
        splited_pole1(5)=splited_pole1(5) &
             +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count)) &
             *(pow1_c1**2) /(alpha_psi1**2) &
             +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
             *(pow1_c0*pow1_c1) /(alpha_psi1**2) &
             +xaxis(i_count)*(1*alpha_psi1**2-charge_density_nume1(i_count))&
             *(pow1_c0**2) /(alpha_psi1**2)
     end do

     if( POWER==2 ) then
        do i_count=1,max_number_xi
           splited_sz2(5,i_count)=splited_sz2(5,i_count) &
                +( sz_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*sz_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*sz_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2) &
                +sz_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*sz_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +sz_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )

           splited_charge2(5,i_count)=splited_charge2(5,i_count) &
                +( charge_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2 ) &
                +charge_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2 ) &
                +2*charge_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2 ) &
                +charge_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
        end do

        do i_count=1,TOTAL_SITE_NUMBER
           splited_charge_density2(5,i_count)= &
                splited_charge_density2(5,i_count) &
                +( charge_density_nume5(i_count)*(pow2_c2**2)/(alpha_psi1**2) &
                +2*charge_density_nume4(i_count)*(pow2_c1*pow2_c2)/(alpha_psi1**2) &
                +2*charge_density_nume3_2(i_count)*(pow2_c0*pow2_c2)/(alpha_psi1**2)&
                +charge_density_nume3_1(i_count)*(pow2_c1**2)/(alpha_psi1**2) &
                +2*charge_density_nume2(i_count)*(pow2_c0*pow2_c1)/(alpha_psi1**2) &
                +charge_density_nume1(i_count)*(pow2_c0**2)/(alpha_psi1**2) )
           splited_pole2(5)=splited_pole2(5) &
                +xaxis(i_count)*(1*alpha_h_h_psi1**2-charge_density_nume5(i_count))&
                *(pow2_c2**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_h_psi1*alpha_h_h_psi1-charge_density_nume4(i_count))&
                *(pow2_c1*pow2_c2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_h_psi1-charge_density_nume3_2(i_count))&
                *(pow2_c0*pow2_c2) /(alpha_psi1**2)&
                +xaxis(i_count)*(1*alpha_h_psi1*alpha_h_psi1-charge_density_nume3_1(i_count))&
                *(pow2_c1**2) /(alpha_psi1**2) &
                +2*xaxis(i_count)*(1*alpha_psi1*alpha_h_psi1-charge_density_nume2(i_count))&
                *(pow2_c0*pow2_c1) /(alpha_psi1**2) &
                +xaxis(i_count)*(1*alpha_psi1*alpha_psi1-charge_density_nume1(i_count))&
                *(pow2_c0**2) /(alpha_psi1**2) 
        end do
     end if
  end if

end subroutine input_splited_cor_number_power
