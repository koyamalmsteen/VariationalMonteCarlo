!
! $Id: sxsy.f90,v 1.2 2003/11/05 09:16:28 k-yukino Exp $
!
#include "parameter.h"

subroutine sxsy(input,tmp_lambda,left_vector_up,left_vector_down,&
     aaa,bbb,tmp_result_sxsy,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: input,tmp_lambda
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: left_vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: left_vector_down
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: tmp_result_sxsy
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer :: i_count,k_count,tmp_up,tmp_down
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_store,global_site_table_down_store 
  real(8),external :: jastrow
  real(8) :: q_up,q_down
  real(8) :: temp_a,temp_b
  integer :: tmp
! init
  tmp_result_sxsy=0
  name="sxsy" ; error=0
! init local
  q_up=0 ; q_down=0
  tmp_up=-1 ; tmp_down=-1
  temp_a=0 ; temp_b=0
  tmp=0

  global_site_table_up_store=global_site_table_up
  global_site_table_down_store=global_site_table_down

  if( aaa==bbb ) then                            ! ナンバーオペレーター
! 第1項
     if( global_site_table_up_store(aaa)==1 ) then
        temp_a=input*tmp_lambda
        if( global_site_table_down_store(aaa)==1 ) then
           temp_a=0
        end if
     end if

! 第2項
     if( global_site_table_down_store(aaa)==1 ) then
        temp_b=input*tmp_lambda
        if( global_site_table_up_store(aaa)==1 ) then
           temp_b=0
        end if
     end if

  else                                        ! aaa/=bbb

! 第1項の計算 (aaa;i_count bbb;j_count)
     if( global_site_table_up_store(aaa)==0 .and. &
          global_site_table_up_store(bbb)==1 .and. &
          global_site_table_down_store(bbb)==0 .and. &
          global_site_table_down_store(aaa)==1 ) then

        gamma_up=left_vector_up
        gamma_down=left_vector_down

! ↑置き換え作業
        do i_count=1,TOTAL_UP_ELECTRON
           if( gamma_up(i_count)==bbb ) then
              gamma_up(i_count)=aaa
              tmp=i_count
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
        do i_count=1,TOTAL_DOWN_ELECTRON
           if( gamma_down(i_count)==aaa ) then
              gamma_down(i_count)=bbb
              tmp=i_count
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

        temp_a=-q_up*q_down*input*jastrow(gamma_up,gamma_down)
     else
        temp_a=0
     end if

! 第2項の計算 (aaa;i_count bbb;j_count)
     if( global_site_table_up_store(bbb)==0 .and. &
          global_site_table_up_store(aaa)==1 .and. &
          global_site_table_down_store(aaa)==0 .and. &
          global_site_table_down_store(bbb)==1 )then   ! もし移動先が空なら移動させる
        
        gamma_up=left_vector_up
        gamma_down=left_vector_down
! ↑置き換え作業
        do i_count=1,TOTAL_UP_ELECTRON
           if( gamma_up(i_count)==aaa ) then
              gamma_up(i_count)=bbb
              tmp=i_count
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
        do i_count=1,TOTAL_DOWN_ELECTRON
           if( gamma_down(i_count)==bbb ) then
              gamma_down(i_count)=aaa
              tmp=i_count
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

  tmp_result_sxsy=0.5*(temp_a+temp_b)

  global_site_table_up=global_site_table_up_store
  global_site_table_down=global_site_table_down_store

end subroutine sxsy
