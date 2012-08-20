!
! $Id: super.f90,v 1.3 2003/11/05 09:17:06 k-yukino Exp $
!
#include "parameter.h"

!
! 方向=1 : 右
! 方向=2 : 下
! 方向=3 : 左
! 方向=4 : 上
!

subroutine super(input,left_vector_up,left_vector_down,aaa,bbb,&
     tmp_result_super,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: input
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: left_vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: left_vector_down
  integer,intent(in) :: aaa,bbb
  real(8),intent(out) :: tmp_result_super
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(TOTAL_SITE_NUMBER) :: global_site_table_up_store,&
       global_site_table_down_store
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  real(8) :: q_up,q_down
  real(8) :: temp_a,temp_b,temp_c,temp_d
  integer :: tmp
  real(8),external :: jastrow
  integer :: i_count,j_count,k_count
  integer :: fugou
! init
  tmp_result_super=0
  name="super" ; error=0
! init local
  global_site_table_up_store=0 ; global_site_table_down_store=0 
  gamma_up=0 ; gamma_down=0
  q_up=0 ; q_down=0
  temp_a=0 ; temp_b=0 ; temp_c=0 ; temp_d=0 ; tmp=0
  fugou=0

  global_site_table_up_store=global_site_table_up
  global_site_table_down_store=global_site_table_down

  if( JIGEN==1 ) then
     error=-1 ; return
  else if( JIGEN==2 ) then
!
! i_count,j_countを固定する
!
     do i_count=1,4
        do j_count=1,4

           fugou=(-1)**(i_count+j_count)

           temp_a=0 ; temp_b=0 ; temp_c=0 ; temp_d=0
           if( neighbor_table2(aaa,i_count)/=0 .and. &
                neighbor_table2(bbb,j_count)/=0 ) then

!!
!! 第1項の計算 (aaa;i_count bbb;j_count)
!!
              if( aaa==bbb ) then              ! upがナンバーオペレーター
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(aaa)==1 ) then
                    q_up=1
                 else
                    q_up=0
                 end if
              else                             ! 通常の演算子
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(aaa)==0 .and. &
                      global_site_table_up_store(bbb)==1 ) then
! ↑置き換え作業
                    tmp=0
                    
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
                       q_up=q_up+d_tilde_up_inverse(k_count,tmp) &
                            *d_tilde_up(tmp,k_count)
                    end do
                    
                 else
                    q_up=0
                 end if
              end if

              if( neighbor_table2(aaa,i_count)== &         ! downがナンバー
                   neighbor_table2(bbb,j_count) ) then

                 gamma_down=left_vector_down

                 if( global_site_table_down_store&
                      &(neighbor_table2(aaa,i_count))==1 ) then
                    q_down=1
                 else
                    q_down=0
                 end if
              else                                         ! 通常の演算子 
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(neighbor_table2(aaa,i_count)&
                      &)==0 .and. global_site_table_down_store(&
                      &neighbor_table(bbb,j_count))==1 ) then



! ↓置き換え作業
                    tmp=0

                    do k_count=1,TOTAL_DOWN_ELECTRON
                       if( gamma_down(k_count)==&
                            neighbor_table2(bbb,j_count) ) then
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
                 else
                    q_down=0
                 end if
              end if
!
! 第1項の結果
!
              temp_a=q_up*q_down*input*jastrow(gamma_up,gamma_down)
              
!!
!! 第2項の計算
!!
              if( aaa==neighbor_table2(bbb,j_count) ) then ! upがナンバー
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(aaa)==1 ) then
                    q_up=1
                 else
                    q_up=0
                 end if
              else                                     ! 通常の演算子
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(aaa)==0 .and. &
                      global_site_table_up_store&
                      &(neighbor_table2(bbb,j_count))==1 ) then
! ↑置き換え作業
                    tmp=0

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
                 else
                    q_up=0
                 end if
              end if
                  
              if( bbb==neighbor_table2(aaa,i_count) ) then ! downがナンバー
                 gamma_down=left_vector_down

                 if(  global_site_table_down_store(bbb)==1 )then 
                    q_down=1
                 else
                    q_down=0
                 end if
              else                                         ! 通常の演算子
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(bbb)==1 .and. &
                      global_site_table_down_store&
                      &(neighbor_table2(aaa,i_count))==0 ) then
! ↓置き換え作業
                    tmp=0

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
                 else
                    q_down=0
                 end if
              end if
!
! 第2項の結果
!
              temp_b=q_up*q_down*input*jastrow(gamma_up,gamma_down)

!!
!! 第3項の計算 (aaa;i_count bbb;j_count)
!!

              if( aaa==bbb ) then              ! downがナンバーオペレーター
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(aaa)==1 ) then
                    q_down=1
                 else
                    q_down=0
                 end if
              else                             ! 通常の演算子
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(aaa)==0 .and. &
                      global_site_table_down_store(bbb)==1 ) then
! ↓置き換え作業
                    tmp=0

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
                 else
                    q_down=0
                 end if
              end if

              if( neighbor_table2(aaa,i_count)== &         ! upがナンバー
                   neighbor_table2(bbb,j_count) ) then
                 gamma_up=left_vector_up

                 if( global_site_table_up_store&
                      &(neighbor_table2(aaa,i_count))==1 ) then
                    q_up=1
                 else
                    q_up=0
                 end if
              else                                         ! 通常の演算子 
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(neighbor_table2(aaa,i_count)&
                      &)==0 .and. global_site_table_up_store(&
                      &neighbor_table(bbb,j_count))==1 ) then
! ↑置き換え作業
                    tmp=0

                    do k_count=1,TOTAL_UP_ELECTRON
                       if( gamma_up(k_count)==&
                            neighbor_table2(bbb,j_count) ) then
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
                 else
                    q_up=0
                 end if
              end if
!
! 第3項の結果
!
              temp_c=q_up*q_down*input*jastrow(gamma_up,gamma_down)
!!
!! 第4項の計算
!!
              if( aaa==neighbor_table2(bbb,j_count) ) then ! downがナンバー
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(aaa)==1 ) then
                    q_down=1
                 else
                    q_down=0
                 end if
              else                                     ! 通常の演算子
                 gamma_down=left_vector_down

                 if( global_site_table_down_store(aaa)==0 .and. &
                      global_site_table_down_store&
                      &(neighbor_table2(bbb,j_count))==1 ) then
! ↓置き換え作業
                    tmp=0

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
                 else
                    q_down=0
                 end if
              end if
                  
              if( bbb==neighbor_table2(aaa,i_count) ) then ! upがナンバー
                 gamma_up=left_vector_up
                    
                 if(  global_site_table_up_store(bbb)==1 )then 
                    q_up=1
                 else
                    q_up=0
                 end if
              else                                           ! 通常の演算子
                 gamma_up=left_vector_up

                 if( global_site_table_up_store(neighbor_table2(aaa,i_count)&
                      &)==0 .and. global_site_table_up_store(&
                      &neighbor_table(bbb,j_count))==1 ) then
! ↑置き換え作業
                    tmp=0

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
                 else
                    q_up=0
                 end if
              end if
!
! 第4項の結果
!
              temp_d=q_up*q_down*input*jastrow(gamma_up,gamma_down)

              tmp_result_super=tmp_result_super &
                   +(temp_a+temp_b+temp_c+temp_d)*fugou*0.25*0.25
           end if
        end do
     end do
     
  end if

  global_site_table_up=global_site_table_up_store
  global_site_table_down=global_site_table_down_store

end subroutine super
