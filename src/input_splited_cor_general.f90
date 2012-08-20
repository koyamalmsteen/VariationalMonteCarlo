!
! $Id: input_splited_cor_general.f90,v 1.2 2004/03/10 03:04:02 k-yukino Exp $
!
#include "parameter.h"

subroutine input_splited_cor_general(sample_number,alpha_psi1,&
     tmp_sxsy0,tmp_super0,name,error)
  use global_variables
  implicit none

  integer,intent(in) :: sample_number
  real(8),intent(in) :: alpha_psi1 
  real(8),dimension(max_number_xi),intent(in) :: tmp_sxsy0,tmp_super0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
! init
  name="input_splited_cor_general" ; error=0

  if( mod(sample_number,5)==1 ) then
     do i_count=1,max_number_xi
        splited_sxsy0(1,i_count)=splited_sxsy0(1,i_count) &
             +tmp_sxsy0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_super0(1,i_count)=splited_super0(1,i_count) &
            +tmp_super0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==2 ) then
     do i_count=1,max_number_xi
        splited_sxsy0(2,i_count)=splited_sxsy0(2,i_count) & 
             +tmp_sxsy0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_super0(2,i_count)=splited_super0(2,i_count) &
            +tmp_super0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==3 ) then
     do i_count=1,max_number_xi
        splited_sxsy0(3,i_count)=splited_sxsy0(3,i_count) &
             +tmp_sxsy0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_super0(3,i_count)=splited_super0(3,i_count) &
            +tmp_super0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==4 ) then
     do i_count=1,max_number_xi
        splited_sxsy0(4,i_count)=splited_sxsy0(4,i_count) &
             +tmp_sxsy0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_super0(4,i_count)=splited_super0(4,i_count) &
            +tmp_super0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  else if( mod(sample_number,5)==0 ) then
     do i_count=1,max_number_xi
        splited_sxsy0(5,i_count)=splited_sxsy0(5,i_count) &
             +tmp_sxsy0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
        splited_super0(5,i_count)=splited_super0(5,i_count) &
            +tmp_super0(i_count)/alpha_psi1/dble(dble(MAX_MONTECARLO_SAMPLE)/5)
     end do
  end if

end subroutine input_splited_cor_general
