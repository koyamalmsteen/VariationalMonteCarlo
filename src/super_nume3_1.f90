!
! $Id: super_nume3_1.f90,v 1.4 2003/11/20 14:02:56 k-yukino Exp $
!
#include "parameter.h"

subroutine super_nume3_1(alpha_h_psi1,vector_up,vector_down,aaa,bbb,result23,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: alpha_h_psi1
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: result23
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: q_up,q_down
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_store,&
       global_site_table_down_store 
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_in&
       ,global_site_table_down_in
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer :: i_count,j_count,k_count
  real(8),external :: jastrow
  real(8) :: temp_a,temp_b,temp_c,temp_d
  real(8) :: result2_1,result2_2,result2_3,result2_4
  real(8) :: result23_1,result23_2
  real(8) :: tmp_result23
  integer :: tmp
  real(8),dimension(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON) :: d_tilde_up_inverse_store
  real(8),dimension(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON) :: d_tilde_down_inverse_store
  real(8) :: input
! init
  name="super_nume3_1" ; error=0
! init local
  temp_a=0 ; temp_b=0
  tmp=0
  result23=0 ; result23_1=0 ; result23_2=0

  global_site_table_up_in=global_site_table_up
  global_site_table_down_in=global_site_table_down
  d_tilde_up_inverse_store=d_tilde_up_inverse
  d_tilde_down_inverse_store=d_tilde_down_inverse

!
! <psi|H|alpha><alpha|O|beta_1><beta_1|H|psi>
!               ^^^^^^^^^^^^^   ^^^^^^^^^^^^^^
!                 result2         result3     
!

!
! result23を求める
!
! <alpha|と結び付く|beta_1>で計算を行なう。
!
  do i_count=1,4
     do j_count=1,4
        if( neighbor_table2(aaa,i_count)/=0 .and. &
             neighbor_table2(bbb,j_count)/=0 ) then

!!
!! 第1項
!!

! ナンバーオペレーターの場合
           if( aaa==bbb .and. &
                neighbor_table2(aaa,i_count)== &
                neighbor_table2(bbb,j_count) ) then

              if( global_site_table_up_store(aaa)==1 .and.&
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==1 &
                   ) then
                 result2_1=1
              end if
           else
! ナンバーオペレーター以外
              if( global_site_table_up_store(aaa)==0 .and.&
                   global_site_table_up_store(bbb)==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(bbb,j_count))==1 &
                   ) then ! もし移動先が空なら移動させる        

                 gamma_up=vector_up
                 gamma_down=vector_down
!
! ↓置き換え作業
!
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==neighbor_table2(bbb,j_count) ) then
                       gamma_down(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,2,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_down=0

                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,tmp) &
                         *d_tilde_down(tmp,k_count)
                 end do
!! 
!! ↑置き換え作業
!!
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==bbb ) then
                       gamma_up(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_up=0

                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_down+d_tilde_up_inverse(k_count,tmp) &
                         *d_tilde_up(tmp,k_count)
                 end do

                 temp_c=-q_up*q_down*input*jastrow(gamma_up,gamma_down)
              else
                 temp_c=0
              end if
           end if

!! 第2項の計算

! ナンバーオペレーターの場合
           if( aaa==neighbor_table2(bbb,j_count) .and. &
                bbb==neighbor_Table2(aaa,i_count) ) then
              
              if( global_site_table_up_store(aaa)==1 .and. &
                   global_site_table_down_store(bbb)==1 )then 
                 result2_2=1
              end if
           else                  
! ナンバーオペレーター以外
              if( global_site_table_up_store(aaa)==0 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(bbb,j_count))==1 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_down_store(bbb)==1 )then 
! もし移動先が空なら移動させる

                 gamma_up=vector_up
                 gamma_down=vector_down
! ↑置き換え作業
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_up(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_up=0
                 
                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_up+d_tilde_up_inverse(k_count,tmp) &
                         *d_tilde_up(tmp,k_count)
                 end do

! ↓置き換え作業
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==bbb ) then
                       gamma_down(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_down=0

                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,tmp) &
                         *d_tilde_down(tmp,k_count)
                 end do

                 temp_b=-q_up*q_down*input*jastrow(gamma_up,gamma_down)
              else
                 temp_b=0
              end if
           end if
!!
!! 第3項
!!

! ナンバーオペレーターの場合
           if( aaa==bbb .and. &
                neighbor_table2(aaa,i_count)== &
                neighbor_table2(bbb,j_count) ) then

              if( global_site_table_down_store(aaa)==1 .and.&
                   global_site_table_up_store&
                   &(neighbor_table2(aaa,i_count))==1 &
                   ) then
                 result2_3=1
              end if
           else
! ナンバーオペレーター以外
              if( global_site_table_down_store(aaa)==0 .and.&
                   global_site_table_down_store(bbb)==1 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(bbb,j_count))==1 &
                   ) then ! もし移動先が空なら移動させる        

                 gamma_up=vector_up
                 gamma_down=vector_down
!
! ↑置き換え作業
!
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==neighbor_table2(bbb,j_count) ) then
                       gamma_up(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_up=0

                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_up+d_tilde_up_inverse(k_count,tmp) &
                         *d_tilde_up(tmp,k_count)
                 end do
!! 
!! ↓置き換え作業
!!
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==bbb ) then
                       gamma_down(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_down=0

                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,tmp) &
                         *d_tilde_down(tmp,k_count)
                 end do

                 temp_c=-q_up*q_down*input*jastrow(gamma_up,gamma_down)
              else
                 temp_c=0
              end if
           end if
!!
!! 第4項の計算
!!

! ナンバーオペレーターの場合
           if( aaa==neighbor_table2(bbb,j_count) .and. &
                bbb==neighbor_Table2(aaa,i_count) ) then

              if( global_site_table_down_store(aaa)==1 .and. &
                   global_site_table_up_store(bbb)==1 )then 
                 result2_4=1
              end if
           else
! ナンバーオペレーター以外
              if( global_site_table_down_store(aaa)==0 .and. &
                   global_site_table_down_store&
                   &(neighbor_table2(bbb,j_count))==1 .and. &
                   global_site_table_up_store&
                   &(neighbor_table2(aaa,i_count))==0 .and. &
                   global_site_table_up_store(bbb)==1 )then 
                     
                 gamma_up=vector_up
                 gamma_down=vector_down

! ↑置き換え作業
                 do k_count=1,TOTAL_UP_ELECTRON
                    if( gamma_up(k_count)==bbb ) then
                       gamma_up(k_count)=neighbor_table2(aaa,i_count)
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_up=0
                     
                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_up+d_tilde_up_inverse(k_count,tmp) &
                         *d_tilde_up(tmp,k_count)
                 end do

! ↓置き換え作業
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    if( gamma_down(k_count)==&
                         neighbor_table2(bbb,j_count) ) then
                       gamma_down(k_count)=aaa
                       tmp=k_count
                       exit
                    end if
                 end do

                 call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
                 call error_check(name,error)
!
! 要素計算と内積
!
                 q_down=0
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,tmp) &
                         *d_tilde_down(tmp,k_count)
                 end do

                 temp_d=-q_up*q_down*input*jastrow(gamma_up,gamma_down)
              else
                 temp_d=0
              end if
           end if

           tmp_result23=result23+(-1)**(i_count+j_count)&
                *(temp_a+temp_b+temp_c+temp_d)
        end if
     end do
  end do

  result23=0.5*result23

  global_site_table_up=global_site_table_up_in
  global_site_table_down=global_site_table_down_in
  d_tilde_up_inverse=d_tilde_up_inverse_store
  d_tilde_down_inverse=d_tilde_down_inverse_store

end subroutine super_nume3_1
