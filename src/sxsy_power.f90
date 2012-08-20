!
! $Id: sxsy_power.f90,v 1.4 2003/11/20 14:05:22 k-yukino Exp $
!

#include "parameter.h"

subroutine sxsy_power(alpha_psi0,alpha_psi1,alpha_h_psi1,alpha_h_h_psi1,&
     lambda,vector_up,vector_down,nume1,nume2,nume3_1,nume3_2,nume4,nume5,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_psi0,lambda  
  real(8),intent(in) :: alpha_h_h_psi1,alpha_h_psi1,alpha_psi1 
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
  real(8),dimension(max_number_xi),intent(out) :: nume1,nume2,nume3_1,nume3_2,nume4,nume5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: tmp_result1,tmp_result2,tmp_result3_1,tmp_result3_2,tmp_result4,tmp_result5,tmp_result
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_store,global_site_table_down_store
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up_inverse_store
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down_inverse_store

  integer :: i_count,j_count
  real(8),dimension(max_number_xi) :: result
  real(8) :: result23,result234

! init
  name="sxsy_power" 
  error=0
! init
  nume1=0 ; nume2=0 ; nume3_1=0 ; nume3_2=0 ; nume4=0 ; nume5=0
  result=0
! local init
  tmp_result1=0 ; tmp_result2=0 ; tmp_result3_1=0
  result23=0 ; result234=0

  global_site_table_up_store=global_site_table_up
  global_site_table_down_store=global_site_table_down
  d_tilde_up_inverse_store=d_tilde_up_inverse
  d_tilde_down_inverse_store=d_tilde_down_inverse

!
! メインルーチン
!
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
!! <sxsy> 
        call sxsy(alpha_psi0,lambda,vector_up,vector_down,&
             i_count,j_count,tmp_result,name,error)
        call error_check(name,error)
! <sxsy>
        tmp_result1=tmp_result*alpha_psi1

        nume1(distance_sequence(i_count,j_count))= &
             nume1(distance_sequence(i_count,j_count)) &
             +tmp_result1

!! <H sxsy>            (<sxsy>の結果を使う)
        tmp_result2=tmp_result*alpha_h_psi1

        nume2(distance_sequence(i_count,j_count))= &
             nume2(distance_sequence(i_count,j_count)) &
             +tmp_result2

!! <H sxsy H>
        call sxsy_nume3_1(alpha_h_psi1,vector_up,vector_down,i_count,j_count,&
             result23,name,error)
        call error_check(name,error);

        tmp_result3_1=alpha_h_psi1*result23

        nume3_1(distance_sequence(i_count,j_count))= &
             nume3_1(distance_sequence(i_count,j_count)) &
             +tmp_result3_1

        if( POWER==2 ) then
!! <O H^2>   ---> <psi|O|alpha><alpha|H^2|psi>とする (<sxsy>の結果を使う)
           tmp_result3_2=tmp_result*alpha_h_h_psi1

           nume3_2(distance_sequence(i_count,j_count))= &
                nume3_2(distance_sequence(i_count,j_count)) &
                +tmp_result3_2
          
!! <H O H^2> ---> <psi|HO|alpha><alpha|H^2|psi> (result23を使う)
! ややこしい
           tmp_result4=result23*alpha_h_h_psi1

           nume4(distance_sequence(i_count,j_count))= &
                nume4(distance_sequence(i_count,j_count)) &
                +tmp_result4

!! <H^2 O H^2> ---> <psi|H^2|alpha><alpha|OH^2|psi>
           call sxsy_nume5(alpha_h_h_psi1,vector_up,vector_down,&
                i_count,j_count,result234,name,error)
           call error_check(name,error)

           tmp_result5=result234*alpha_h_h_psi1

           nume5(distance_sequence(i_count,j_count))= &
                nume5(distance_sequence(i_count,j_count)) &
                +tmp_result5
        end if

     end do
  end do

  global_site_table_up=global_site_table_up_store
  global_site_table_down=global_site_table_down_store
  d_tilde_up_inverse=d_tilde_up_inverse_store
  d_tilde_down_inverse=d_tilde_down_inverse_store

end subroutine sxsy_power
