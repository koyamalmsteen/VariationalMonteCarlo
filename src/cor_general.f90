!
! $Id: cor_general.f90,v 1.3 2003/11/18 08:05:41 k-yukino Exp $
!
#include "parameter.h"

subroutine cor_general(alpha_psi0,lambda,vector_up,vector_down,&
     tmp_sxsy0,tmp_super0,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_psi0,lambda
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),dimension(max_number_xi),intent(out) :: tmp_sxsy0,&
       tmp_super0
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: tmp_result
  integer ::  i_count,j_count,k_count
! init
  tmp_sxsy0=0 ; tmp_super0=0
  name="cor" ; error=0
  tmp_result=0

!
! ★★★要素計算ルーチン★★★
!
 
! sxsy
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=i_count,TOTAL_SITE_NUMBER ! ダブルカウントをしないため
        tmp_result=0

        call sxsy(alpha_psi0,lambda,vector_up,vector_down,&
             i_count,j_count,tmp_result,name,error)
        call error_check(name,error)
           
        tmp_sxsy0(distance_sequence(i_count,j_count))&
             &=tmp_sxsy0(distance_sequence(i_count,j_count))&
             +tmp_result
     end do
  end do

! super
  if( JIGEN==2 ) then
     do i_count=1,TOTAL_SITE_NUMBER
        do j_count=i_count,TOTAL_SITE_NUMBER ! ダブルカウント防止
           tmp_result=0
              
           call super(alpha_psi0,vector_up,vector_down,i_count,j_count,&
                tmp_result,name,error)
           call error_check(name,error)

           tmp_super0(distance_sequence(i_count,j_count))&
                &=tmp_super0(distance_sequence(i_count,j_count)) &
                +tmp_result
        end do
     end do
  end if

end subroutine cor_general
