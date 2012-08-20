!
! $Id: calc_acceptance_ratio.f90,v 1.5 2004/02/15 10:52:35 k-yukino Exp k-yukino $
!
#include "parameter.h"

subroutine calc_acceptance_ratio(gamma_up_work,gamma_down_work,tmp_lambda,&
                                 denomi_lambda,idou_spin,idou_moto,acceptance_ratio,name,error)
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON),intent(in) :: gamma_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(in) :: gamma_down_work
  real(8),intent(in) :: tmp_lambda
  real(8),intent(out) :: denomi_lambda
  integer,intent(in) :: idou_spin,idou_moto
  real(8),intent(out) :: acceptance_ratio
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count
  real(8),external :: jastrow,lambda_field
  real(8) :: q
! init
  acceptance_ratio=0 ; name="calc_acceptance_ratio" ; error=0
! init [local]
  denomi_lambda=0
  q=0

  if( idou_spin==1 ) then         ! ↑が動いたなら
     do i_count=1,TOTAL_UP_ELECTRON
        q=q+d_tilde_up_inverse(i_count,idou_moto)&
             *d_tilde_up(idou_moto,i_count)
     end do

  else if( idou_spin==2 ) then     ! ↓が動いたなら
     do i_count=1,TOTAL_DOWN_ELECTRON
        q=q+d_tilde_down_inverse(i_count,idou_moto) &
             *d_tilde_down(idou_moto,i_count)
     end do
  end if

! 変分パラメータの数が>0なら
  if( NUMBER_XI>0 ) then
     denomi_lambda=jastrow(gamma_up_work,gamma_down_work)
  end if

! 電場の効果があるなら
  if( electric_field/=0 .and. PERIODIC==1 .and. xi_field/=0 ) then
     denomi_lambda=denomi_lambda*lambda_field(gamma_up_work,gamma_down_work)
  end if

  q=q*denomi_lambda/tmp_lambda

  acceptance_ratio=(abs(q))**2

end subroutine calc_acceptance_ratio



#ifdef DEBUG_CALC_ACCEPTANCE_RATIO
program main
  use global_variables
  implicit none

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up,gamma_up_work
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down,gamma_down_work
  character(32) :: name
  integer :: error
!
  real(8) :: ransuu
  integer :: number_accept,number_reject
! local
  real(8) :: denominator_up,denominator_down,denominator
  integer :: i_count,j_count,k_count
  real(8) :: acceptance_ratio
! local 
  real(8) :: temp
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  real(8) :: bunbo,bunbo_up,bunbo_down
  integer :: idou_spin,idou_moto
  integer,dimension(1) :: ic
  character(32) :: filename
  integer :: kari
! NAG関連
  real(8),dimension(:),allocatable :: wk_space
  integer :: ifail 
! MPI関連
  integer :: IERROR
! init
  denominator_up=0 ; denominator_down=0
  ic(1)=1
!
  zero_approx_energy=0 ; hf_iteration_result=0

  number_accept=0 ; number_reject=0
  d_tilde_up_inverse=0 ; d_tilde_down_inverse=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)

! 初期gammaを乱数より決定。最初だけsortをしておく(別にする必要はないが)
  call choice_gamma_sigma(TOTAL_UP_ELECTRON,gamma_up,name,error)
  call error_check(name,error)
  call choice_gamma_sigma(TOTAL_DOWN_ELECTRON,gamma_down,name,error)
  call error_check(name,error)

!
! さて、この時点で最初の電子状態は決定された。
! (まだサンプル数をインクリメントしません。後でするので。)
!

! さて、 初期の電子状態を調べる。

  k_count=0

  do i_count=1,TOTAL_UP_ELECTRON
     do j_count=1,TOTAL_DOWN_ELECTRON
        if( gamma_up(i_count)==gamma_down(j_count) ) then ! ダブルオキュパイ
           k_count=k_count+1
        end if
     end do
  end do

! 初期γの逆行列を作る
  call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
  call error_check(name,error)
  call make_d_tilde(gamma_up,gamma_down,2,1,name,error)
  call error_check(name,error)
  call invert(1,1,name,error) 
  call error_check(name,error)
  call invert(2,1,name,error) 
  call error_check(name,error)

  acceptance_ratio=1

! 試行電子状態を作業領域に

  gamma_up_work=gamma_up ; gamma_down_work=gamma_down

! ここからメインのループ

  temp=0

!  open(21,file="aho")

  do while( number_accept<MAX_MONTECARLO_SAMPLE )
!
!     if( mod(number_accept,10000)==0 .and. temp/=number_accept ) then
!        write(21,*) number_accept,&
!                  &dble(100*number_accept)/dble(number_accept+number_reject)
!        temp=number_accept
!     end if

     call fortran_random(ransuu,name,error)
     call error_check(name,error)

     if( ransuu<=acceptance_ratio ) then

        number_accept=number_accept+1
        gamma_up=gamma_up_work ; gamma_down=gamma_down_work ! 正式に採用した
!
! ここで、エネルギー期待値計算の分母の<alpha|psi>をしておく。
!        
        bunbo_up=1
        do i_count=1,TOTAL_UP_ELECTRON
           bunbo_up=bunbo_up*d_tilde_up(i_count,i_count)
        end do

        bunbo_down=1
        do i_count=1,TOTAL_DOWN_ELECTRON
           bunbo_down=bunbo_down*d_tilde_down(i_count,i_count)
        end do

        bunbo = bunbo_up * bunbo_down

! 要素計算!!!

        if( idou_spin==1 )then        ! ↑スピンが動いたら
           call make_d_tilde(gamma_up,gamma_down,1,1,name,error)
           call error_check(name,error)

           call invert(1,1,name,error) 
           call error_check(name,error)
        end if
     else                                  ! リジェクトされたら 
        number_reject=number_reject+1
     end if

! 新しい電子状態γ'を作る。
! (古い電子状態はgamma_up,gamma_downとして残しておく(リジェクト時に使うから))

     gamma_up_work=gamma_up ; gamma_down_work=gamma_down ! ワーク用に移す

     call choice_new_gamma(gamma_up_work,gamma_down_work,&
                           idou_spin,&
                           idou_moto,name,error)
     call error_check(name,error)

! 以下で、新しい試行電子状態の部分行列を作る。
! 部分行列を作る処理は移動したスピンのみで行なえばよい。
! 移動してなければ変更はないから。

     if( idou_spin==1 ) then  ! ↑スピンが動いたら
        call make_d_tilde(gamma_up_work,gamma_down_work,1,1,name,error)
        call error_check(name,error)
! これは逆行列を作る必要がない。逆行列と掛け合わせる行列だから
     else if( idou_spin==2 ) then ! ↓スピンが動いたら 
        call make_d_tilde(gamma_up_work,gamma_down_work,2,1,name,error)
        call error_check(name,error)
! これは逆行列を作る必要がない。逆行列と掛け合わせる行列だから
     end if

! アクセプタンス・レシオを求める。

     call calc_acceptance_ratio(gamma_up,gamma_down,gamma_up_work,&
                                gamma_down_work,idou_spin,&
                                idou_moto,acceptance_ratio,name,error)
     call error_check(name,error)
  end do

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CALC_ACCEPTANCE_RATIO */
