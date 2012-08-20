!
! $Id: calc_element_pow5.f90,v 1.4 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!

!
! E=<phi_0 |H^5|phi_0 >
!

#include "parameter.h"

subroutine calc_element_pow5(vector_up,vector_down,alpha_h_h_h_psi1,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(out) :: alpha_h_h_h_psi1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: double_occupy
  real(8) :: input,input_up,input_down
  integer :: sign
  integer,dimension(TOTAL_SITE_NUMBER) :: right_site_table4
  real(8),external :: jastrow,lambda_field
  integer :: i_count,b4u_count,b4d_count,b5u_count,b5d_count
  integer,dimension(TOTAL_UP_ELECTRON) :: right_vector_up4,right_vector_up5
  integer,dimension(TOTAL_DOWN_ELECTRON) :: right_vector_down4,&
                                            right_vector_down5
  real(8) :: result4,result5,result6_7
  integer :: number_ket_vector_up,number_ket_vector_down,&
             number_ket_vector_up5,number_ket_vector_down5
  integer,dimension(2*JIGEN*TOTAL_UP_ELECTRON+1,TOTAL_UP_ELECTRON) :: available_ket_up_table,available_ket_up_table5
  integer,dimension(2*JIGEN*TOTAL_DOWN_ELECTRON+1,TOTAL_DOWN_ELECTRON) :: available_ket_down_table,available_ket_down_table5
  real(8) :: tmp_jastrow
  real(8),external :: calc_total_pole
  real(8) :: pole
! init
  name="calc_element_pow3" ; error=0
! init [local]
  double_occupy=0 
  right_site_table4=0
  right_vector_up4=0 ; right_vector_down4=0 
  right_vector_up5=0 ; right_vector_down5=0
  alpha_h_h_h_psi1=0 ; result4=0 ; result5=0
  number_ket_vector_up=0 ; number_ket_vector_down=0
  number_ket_vector_up5=0 ; number_ket_vector_down5=0
  available_ket_up_table=0 ; available_ket_down_table=0
  available_ket_up_table5=0 ; available_ket_down_table5=0
  input=0 ; input_up=0 ; input_down=0 ; result6_7=0 ; sign=0
  tmp_jastrow=0
  pole=0

!
! <phi_0|H^5|phi_0> --> 
! <phi_0 |beta_1><beta_1|H|beta_2><beta_2|H|alpha><alpha|H|beta_3>
!
! ^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^
!                                                     result4      
!
! <beta_3|H|beta_4><beta_4|H|beta_5><beta_5|phi_0>
!
! ^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^
!  result5              result6        result7
!

! alpha_h_h_h_psi1を求める。
! <alpha|と結び付く|beta_3>を求める。

  call choice_ket_vector(vector_up,vector_down,&
                         number_ket_vector_up,number_ket_vector_down,&
                         available_ket_up_table,available_ket_down_table,&
                         name,error)
  call error_check(name,error)

! result4
  do b4u_count=1,number_ket_vector_up

           do i_count=1,TOTAL_UP_ELECTRON
              right_vector_up4(i_count)=available_ket_up_table(b4u_count,&
                   i_count)
           end do

     do b4d_count=1,number_ket_vector_down
! 電子を一つしか動かせないのでどちらかが必ず固定される
        if( b4u_count==number_ket_vector_up .or. &
            b4d_count==number_ket_vector_down ) then

/*
           do i_count=1,TOTAL_UP_ELECTRON
              right_vector_up4(i_count)=available_ket_up_table(b4u_count,&
                   i_count)
           end do
*/
           do i_count=1,TOTAL_DOWN_ELECTRON
              right_vector_down4(i_count)=available_ket_down_table(&
                   b4d_count,i_count)
           end do

! ケットと全く同じものではさんだときクーロンが作用

           if( b4u_count==number_ket_vector_up .and. &   
               b4d_count==number_ket_vector_down ) then

! ダブルオキュパイの数を調べる。
! その為にサイト表示にする。

              right_site_table4=0
              do i_count=1,TOTAL_UP_ELECTRON
                 right_site_table4(right_vector_up4(i_count))=1
              end do
              do i_count=1,TOTAL_DOWN_ELECTRON
                 right_site_table4(right_vector_down4(i_count))=&
                 right_site_table4(right_vector_down4(i_count))+1
              end do

              double_occupy=0

              do i_count=1,TOTAL_SITE_NUMBER
                 if( right_site_table4(i_count)==2 ) then
                    double_occupy=double_occupy+1
                 end if
              end do

              result4=double_occupy*dble(COULOMB)

!
! 電界の項
!
              if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
                 pole=calc_total_pole(right_vector_up4,right_vector_down4)

                 result4=result4+electric_field*pole
              end if
           else
              result4=dble(TRANSFER)
           end if
!
! result5
!
           call choice_ket_vector(right_vector_up4,right_vector_down4,&
                number_ket_vector_up5,number_ket_vector_down5,&
                available_ket_up_table5,available_ket_down_table5,&
                name,error)
           call error_check(name,error)

           do b5u_count=1,number_ket_vector_up5

              do i_count=1,TOTAL_UP_ELECTRON
                 right_vector_up5(i_count)= &
                      available_ket_up_table5(b5u_count,i_count)
              end do

! jastrow関数の計算を高速化させるため、以下のようにサイト表示に書き換える
! (書き換え作業を考えても高速化されると思う)(calc_elementで使う。必須)
              global_site_table_up=0
              do i_count=1,TOTAL_UP_ELECTRON
                 global_site_table_up(right_vector_up5(i_count))=1
              end do
! result6_7
! <right_vector_up5|と確実に結び付く|right_vector_up5>の逆行列を計算しておく(高速化)

              call make_d_tilde(right_vector_up5,right_vector_down5,1,2,&
                   name,error)
              call error_check(name,error)
              call invert(1,2,name,error)
              call error_check(name,error)
!いんなー
              input_up=1
              sign=0
              
              do i_count=1,TOTAL_UP_ELECTRON
                 input_up=input_up*after_lu_up(i_count,i_count)
                 if( ipiv_up(i_count)/=i_count ) then
                    sign=sign+1
                 end if
              end do

              input_up=input_up*((-1)**sign)

              do b5d_count=1,number_ket_vector_down5
                 
                 if( b5u_count==number_ket_vector_up5 .or. &  ! 電子を一つしか 
                      b5d_count==number_ket_vector_down5 ) then ! 動かせない

                    do i_count=1,TOTAL_DOWN_ELECTRON
                       right_vector_down5(i_count)= &
                            available_ket_down_table5(b5d_count,i_count)
                    end do

! jastrow関数の計算を高速化させるため、以下のようにサイト表示に書き換える
! (書き換え作業を考えても高速化されると思う)(calc_elementで使う。必須)
                    global_site_table_down=0
                    do i_count=1,TOTAL_DOWN_ELECTRON
                       global_site_table_down(right_vector_down5(i_count))=1
                    end do

! <内積計算> (down分)
                    call make_d_tilde(right_vector_up5,right_vector_down5,&
                                      2,2,name,error)
                    call error_check(name,error)
                    call invert(2,2,name,error)
                    call error_check(name,error)

                    input_down=1
                    sign=0
                    do i_count=1,TOTAL_DOWN_ELECTRON
                       input_down=input_down*after_lu_down(i_count,i_count)
                       if( ipiv_down(i_count)/=i_count ) then
                          sign=sign+1
                       end if
                    end do
                    
                    input_down=input_down*((-1)**sign)


                    if( b5u_count==number_ket_vector_up5 .and. &
                         b5d_count==number_ket_vector_down5 ) then

! ダブルオキュパイの数を調べる。
                       double_occupy=0

                       do i_count=1,TOTAL_SITE_NUMBER
                          if( global_site_table_up(i_count)==1 .and. &
                               global_site_table_down(i_count)==1 ) then 
                             double_occupy=double_occupy+1
                          end if
                       end do

                       result5=double_occupy*dble(COULOMB)
                    else
                       result5=dble(TRANSFER)
                    end if

                    tmp_jastrow=jastrow(right_vector_up5,right_vector_down5)

                    if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
                       tmp_jastrow=tmp_jastrow*lambda_field(right_vector_up5,right_vector_down5)
                    end if

                    input = input_up * input_down

                    call calc_element(input,tmp_jastrow,right_vector_up5,&
                         right_vector_down5,result6_7,&
                         name,error)
                    call error_check(name,error)

                    alpha_h_h_h_psi1=alpha_h_h_h_psi1+result4*result5*result6_7

                 end if
              end do
           end do
        end if

     end do
  end do

end subroutine calc_element_pow5
