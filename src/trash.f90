!
! $Id: trash.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine trash(vector_up,vector_down,name,error)
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: vector_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: vector_down
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(TOTAL_UP_ELECTRON) :: vector_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON) :: vector_down_work
  integer :: number_accept,idou_spin,idou_saki,idou_moto
  real(8) :: acceptance_ratio,alpha_psi,alpha_psi_up,alpha_psi_down,&
             tmp_lambda,ransuu
  integer :: sign 
  real(8),external :: jastrow
  integer :: i_count
  real(8) :: denomi_lambda
! init
  name="trash" ; error=0
! 
  number_accept=0 ; idou_spin=0 ; idou_saki=0
  idou_moto=0 ; tmp_lambda=0 ; sign=0 ; ransuu=0
  alpha_psi=0 ; alpha_psi_up=0 ; alpha_psi_down=0
  denomi_lambda=0

  acceptance_ratio=1         ! 最初の状態は必ずアクセプトするので

! TRASHモンテカルロステップだけ電子状態を捨てる。

  vector_up_work=vector_up ; vector_down_work=vector_down
  number_accept=1

  do while( number_accept<=TRASH*(TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON) )
! (ここでモンテカルロサンプルの確率計算をする)
     call fortran_random(ransuu,name,error)
     call error_check(name,error)

     if( ransuu<=acceptance_ratio ) then
        number_accept=number_accept+1
        vector_up=vector_up_work          ! 正式採用の場合更新する
        vector_down=vector_down_work      !
     else
        number_accept=number_accept+1
!!
!!      この場合電子状態に変更はない。前の電子状態をそのまま採り入れる 
!!
     end if

! 高速化のためにサイト表示にしておく
     if( idou_spin==1 ) then
        global_site_table_up=0
        do i_count=1,TOTAL_UP_ELECTRON
           global_site_table_up(vector_up(i_count))=1
        end do
     else
        global_site_table_down=0
        do i_count=1,TOTAL_DOWN_ELECTRON
           global_site_table_down(vector_down(i_count))=1
        end do
     end if

! 逆行列を作る方法
     if( idou_spin==1 )then   ! ↑スピンが動いたら
        call make_d_tilde(vector_up,vector_down,1,1,name,error)
        call error_check(name,error)
        call invert(1,1,name,error) 
        call error_check(name,error)
     else                                   ! ↓スピンが動いたら
        call make_d_tilde(vector_up,vector_down,2,1,name,error)
        call error_check(name,error)
        call invert(2,1,name,error)
        call error_check(name,error)
     end if
!
! ここで、エネルギー期待値計算の分母の<alpha|psi>=lambda_alpha<alpha|HF>
! をしておく。
!
     alpha_psi_up=1 
     sign=0        
     do i_count=1,TOTAL_UP_ELECTRON
        alpha_psi_up=alpha_psi_up*after_lu_up(i_count,i_count)
        if( ipiv_up(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     alpha_psi_up=alpha_psi_up*((-1)**sign)

     alpha_psi_down=1
     sign=0
     do i_count=1,TOTAL_DOWN_ELECTRON
        alpha_psi_down=alpha_psi_down*after_lu_down(i_count,i_count)
        if( ipiv_down(i_count)/=i_count ) then
           sign=sign+1
        end if
     end do

     alpha_psi_down=alpha_psi_down*((-1)**sign)

     tmp_lambda=jastrow(vector_up,vector_down)
     alpha_psi = alpha_psi_up * alpha_psi_down*tmp_lambda

! 新しい電子状態γ'を作る。
! (古い電子状態はvector_up,vector_downとして残しておく(リジェクト時に使うから))

     vector_up_work=vector_up ; vector_down_work=vector_down   ! ワーク用に移す

     if( PROJECTION==0 ) then 
        call choice_new_gamma(vector_up_work,vector_down_work,&
             idou_spin,idou_moto,idou_saki,name,error)
        call error_check(name,error)
     else
        call choice_new_gamma_projection(vector_up_work,vector_down_work,&
             idou_spin,idou_moto,idou_saki,name,error)
        call error_check(name,error)
     end if

! 以下で、新しい試行電子状態の部分行列を作る。
! 部分行列を作る処理は移動したスピンのみで行なえばよい。
! 移動してなければ変更はないから。

! 逆行列を作る方法  
     if( idou_spin==1 ) then  ! ↑スピンが動いたら
        call make_d_tilde(vector_up_work,vector_down_work,1,1,name,error)
        call error_check(name,error)
! downは移動していないので、何もしない
     else if( idou_spin==2 ) then ! ↓スピンが動いたら 
        call make_d_tilde(vector_up_work,vector_down_work,2,1,name,error)
        call error_check(name,error)
! upは移動していないので、何もしない
     end if

! アクセプタンス・レシオの計算
     call calc_acceptance_ratio(vector_up_work,vector_down_work,tmp_lambda,&
                                denomi_lambda,idou_spin,idou_moto,&
                                acceptance_ratio,name,error)
     call error_check(name,error)
  end do

end subroutine trash



#ifdef DEBUG_TRASH
program main
  implici nonte

end program main
#endif /* DEBUG_TRASH */
