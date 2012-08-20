!
! $Id: calc_element.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine calc_element(input,tmp_jastrow,vector_up,vector_down,result,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: input,tmp_jastrow
  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: vector_down
  real(8),intent(out) :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8) :: q_up,q_down
  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  integer :: double_occupy
  real(8),external :: jastrow,calc_total_pole
  real(8) :: pole 
  integer :: i_count,j_count,k_count
  integer :: sign
  real(8) :: result_up,result_down,result_coulomb,result_field
  real,external :: lambda_field
! init
  result=0 ; name="calc_element" ; error=0
! init [local]
  gamma_up=0 ; gamma_down=0 ; q_up=0 ; q_down=0
  double_occupy=0 ; sign=0
  result_up=0 ; result_down=0 ; result_coulomb=0 ; result_field=0
  pole=0

!!
!! トランスファーの↑スピンの項の計算
!!
  gamma_down=vector_down

  do j_count=1,TOTAL_SITE_NUMBER ! 移動先
     do i_count=1,TOTAL_UP_ELECTRON ! 移動元
        if( FLAT_BAND/=3 ) then              ! TASAKI模型以外
           if( neighbor_table(vector_up(i_count),j_count)==1 ) then! ←隣接なら移動可

              if( global_site_table_up(j_count)==0 ) then ! もし移動先が空なら移動させる
                 gamma_up=vector_up
                 gamma_up(i_count)=j_count                ! ここだけ置き換える

                 call make_d_tilde(gamma_up,gamma_down,1,2,name,error) ! workに
                 call error_check(name,error)
! 要素計算と内積
                 q_up=0

                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_up+d_tilde_up_inverse(k_count,i_count) &
                         *d_tilde_up_work(i_count,k_count)
                 end do

                 if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
                    result_up=result_up+dble(TRANSFER)*q_up*input&
                         *jastrow(gamma_up,gamma_down)*lambda_field(gamma_up,gamma_down)
                 else
                    result_up=result_up+dble(TRANSFER)*q_up*input&
                         *jastrow(gamma_up,gamma_down)
                 end if
              end if
           end if
        else                                    ! TASAKI模型の場合
           if( neighbor_table(vector_up(i_count),j_count)==2 ) then
              if( global_site_table_up(j_count)==0 ) then ! もし移動先が空なら移動させる
                 gamma_up=vector_up
                 gamma_up(i_count)=j_count           ! ここだけ置き換える
                 call make_d_tilde(gamma_up,gamma_down,1,2,name,error) ! workに
                 call error_check(name,error)

! 要素計算と内積

                 q_up=0

                 do k_count=1,TOTAL_UP_ELECTRON
                    q_up=q_up+d_tilde_up_inverse(k_count,i_count) &
                         *d_tilde_up_work(i_count,k_count)
                 end do
              
! TASAKI模型の斜めのトランスファーはt/2
                 result_up=result_up+0.5*dble(TRANSFER) &
                      *q_up*input*tmp_jastrow &
                      *jastrow(gamma_up,gamma_down)
              end if
           end if
        end if

     end do
  end do


!!
!! トランスファーの↓の項の計算
!!
  gamma_up=vector_up

  do j_count=1,TOTAL_SITE_NUMBER ! 移動先
     do i_count=1,TOTAL_DOWN_ELECTRON ! 移動元
        if( FLAT_BAND/=3 ) then                  ! TASAKI模型以外
           if( neighbor_table(vector_down(i_count),j_count)==1 ) then! ←隣接なら移動可
              if( global_site_table_down(j_count)==0 ) then ! もし移動先が空なら要素計算と内積

! 要素計算と内積

                 gamma_down=vector_down
                 gamma_down(i_count)=j_count                ! ここだけ置き換える
                 call make_d_tilde(gamma_up,gamma_down,2,2,name,error) ! workに
                 call error_check(name,error)

                 q_down=0

                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,i_count) &
                         *d_tilde_down_work(i_count,k_count)
                 end do

                 if( electric_field/=0 .and. PERIODIC==1 .and.xi_field/=0 ) then
                    result_down=result_down+dble(TRANSFER)*q_down*input&
                         *jastrow(gamma_up,gamma_down) &
                         *lambda_field(gamma_up,gamma_down)
                 else
                    result_down=result_down+dble(TRANSFER)*q_down*input&
                         *jastrow(gamma_up,gamma_down)
                 end if
              end if
           end if
        else                                     ! TASAKI模型の場合
           if( neighbor_table(vector_down(i_count),j_count)==2 ) then
              if( global_site_table_down(j_count)==0 ) then ! もし移動先が空なら要素計算と内積

! 要素計算と内積

                 gamma_down=vector_down
                 gamma_down(i_count)=j_count             ! ここだけ置き換える

                 call make_d_tilde(gamma_up,gamma_down,2,2,name,error) ! workに
                 call error_check(name,error)

                 q_down=0
                 do k_count=1,TOTAL_DOWN_ELECTRON
                    q_down=q_down+d_tilde_down_inverse(k_count,i_count) &
                         *d_tilde_down_work(i_count,k_count)
                 end do

! TASAKI模型の斜めのトランスファーはt/2
                 result_down=result_down+0.5*dble(TRANSFER) &
                      *q_down*input*jastrow(gamma_up,gamma_down)
              end if
           end if
        end if

     end do
  end do

!
! クーロン斥力の項 ( vector_up,down自身ではさんだ時のみ生き残る )
!

! まずダブルオキュパイの数を調べる。
  double_occupy=0

  do i_count=1,TOTAL_SITE_NUMBER
     if( global_site_table_up(i_count)==1 .and. &
         global_site_table_down(i_count)==1 ) then
        double_occupy=double_occupy+1
     end if
  end do

  result_coulomb=dble(COULOMB)*double_occupy*input*tmp_jastrow

!
! vector_up,down自身で挟んだときのみ残る
! 
! 電界の項( EP ) Eはスカラーかつ正
!

  if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
! まずPを求める
     pole=calc_total_pole(vector_up,vector_down)

     result_field=electric_field*pole*input*tmp_jastrow
  end if

  result=result_up+result_down+result_coulomb+result_field

end subroutine calc_element



#ifdef DEBUG_CALC_ELEMENT
program main
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down
  real(8) :: result,tmp_jastrow,input_up,input_down,alpha_psi0
  real(8),external :: jastrow
  character(32) :: name
  integer :: error
! local
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  integer :: i_count
  integer :: sign
  integer,dimension(2) :: ic
! init
  result=0 ; tmp_jastrow=0 ; error=0 ; name="main"
  vector_up=0 ; vector_down=0
  input_up=0 ; input_down=0 ; alpha_psi0=0
  sign=0
  ic=1
  zero_approx_energy=0 ; hf_iteration_result=0

  call MPI_INIT(IERROR)

! パラメーターの初期化
  call init_global_variables(name,error)
  call error_check(name,error)

! 乱数テーブルの初期化
  call random_init(ic,name,error)
  call error_check(name,error)

! 第0近似の実行(INIT_WAVE_VECTOR = 0;Hartree-Fock , 1 ;Plane wave )
! (MYRANK==0のみが実行し、結果をbcastする)
  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "input vector_up" ; read(*,*) vector_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "input vector_down" ; read(*,*) vector_down(i_count)
  end do

  call make_d_tilde(vector_up,vector_down,1,1,name,error)
  call error_check(name,error)
  call make_d_tilde(vector_up,vector_down,2,1,name,error)
  call error_check(name,error)
! 逆行列を求める
  call invert(1,1,name,error)
  call error_check(name,error)
  call invert(2,1,name,error)
  call error_check(name,error)
! <alpha|psi0>を求める

  input_up=1
  sign=0
  do i_count=1,TOTAL_UP_ELECTRON
     input_up=input_up*after_lu_up(i_count,i_count)
     if( ipiv_up(i_count)/=i_count )then
        sign=sign+1
     end if
  end do

  input_up=input_up*((-1)**sign)

  write(*,*) "input_up=",input_up

  input_down=1
  sign=0
  do i_count=1,TOTAL_DOWN_ELECTRON
     input_down=input_down*after_lu_down(i_count,i_count)
     if( ipiv_down(i_count)/=i_count )then
        sign=sign+1
     end if
  end do

  input_down=input_down*((-1)**sign)

  write(*,*) "input_down=",input_down

  alpha_psi0 = input_up * input_down
  tmp_jastrow=jastrow(vector_up,vector_down)

  global_site_table_up=0
  do i_count=1,TOTAL_UP_ELECTRON
     global_site_table_up(vector_up(i_count))=1
  end do

  global_site_table_down=0
  do i_count=1,TOTAL_DOWN_ELECTRON
     global_site_table_down(vector_down(i_count))=1
  end do

  write(*,*) "alpha_psi0=",alpha_psi0
  write(*,*) "tmp_jastrow=",tmp_jastrow
  write(*,*) "vector_up=",vector_up
  write(*,*) "vector_down=",vector_down
  write(*,*) 
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*) 

  call calc_element(alpha_psi0,tmp_jastrow,vector_up,&
                         vector_down,result,name,error)

  write(*,*) "result=",result

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ELEMENT */ 
