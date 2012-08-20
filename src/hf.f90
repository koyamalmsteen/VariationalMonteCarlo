!
! $Id: hf.f90,v 1.26 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine hf(hf_energy,hf_iteration_result,name,error)
  use global_variables
  implicit none

  real(8),intent(out) :: hf_energy
  integer,intent(out) :: hf_iteration_result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8),dimension(:,:),allocatable :: rho_up,rho_down,rho_up_new,rho_down_new
  real(8) :: ransuu,ransuu_tmp,result_up,result_down
  integer :: ier
  integer :: i_count,j_count
  real(8) :: hf_energy_coulomb,hf_energy_transfer
! lapack
  integer :: lapack_nnn,lapack_lda,lapack_lwork,lapack_ifail
  real(8),dimension(:),allocatable :: lapack_www_up,lapack_www_down
  real(8),dimension(:),allocatable :: lapack_work
  character*1 :: lapack_job,lapack_uplo
! init
  hf_energy=0 ; hf_iteration_result=0
  name="hf" ; error=0
! int local
  i_count=0 ; j_count=0 ; result_up=0 ; result_down=0
  hf_energy_transfer=0 ; hf_energy_coulomb=0

! init
  lapack_job='v'                   ! eigen value and eigen vector are needed
  lapack_uplo='u'                  ! upper triangular part of A is stored
  lapack_nnn=TOTAL_SITE_NUMBER
  lapack_lda=TOTAL_SITE_NUMBER
  lapack_lwork=64*TOTAL_SITE_NUMBER
  lapack_ifail=0

  allocate(lapack_www_up(lapack_lda),stat=ier)
  call stat_check("lapack_www_up","hf",1,ier)
  allocate(lapack_www_down(lapack_lda),stat=ier)
  call stat_check("lapack_www_down","hf",1,ier)
  allocate(lapack_work(lapack_lwork),stat=ier)
  call stat_check("lapack_work","hf",1,ier)
! init
  lapack_www_up=0 ; lapack_www_down=0 ; lapack_work=0
!
  allocate(rho_up(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_up","hf",1,ier)
  allocate(rho_down(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_down","hf",1,ier)
  allocate(rho_up_new(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_up_new","hf",1,ier)
  allocate(rho_down_new(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("rho_down_new","hf",1,ier)
! init
  unitary_d_up=0 ; unitary_d_down=0
  rho_up=0 ; rho_down=0
  rho_up_new=0 ; rho_down_new=0

!!! decision of primitive rho_up and rho_down using random number

!
! ☆１ 初期のrhoを決めるルーチン(乱数テーブルの初期化はメインで)
!
  if( HF_INITIAL_RHO==0 ) then                 ! 1,0,1,0...  
     call fortran_random(ransuu_tmp,name,error)! ↑↓か↓↑をを決めるだけ
     call error_check(name,error)              ! 一番低エネルギーが得られる
                                               ! おすすめ！
     do i_count=1,TOTAL_SITE_NUMBER
        if( mod(i_count,2)==1 ) then
           if( ransuu_tmp<0.5 ) then
              rho_up(i_count,i_count)=1 ; rho_down(i_count,i_count)=0 
           else
              rho_up(i_count,i_count)=0 ; rho_down(i_count,i_count)=1
           end if
        else
           if( ransuu_tmp<0.5 ) then
              rho_up(i_count,i_count)=0 ; rho_down(i_count,i_count)=1
           else
              rho_up(i_count,i_count)=1 ; rho_down(i_count,i_count)=0
           end if
        endif
     end do
  else if( HF_INITIAL_RHO==1 ) then     ! sdw的な解を仮定する
                                        ! 0.5>,0.5<,0.5>,0.5<(もしくは逆)と
                                        ! なっている
     call fortran_random(ransuu_tmp,name,error) ! 並べる順番を決める為
     call error_check(name,error)

     do i_count=1,TOTAL_SITE_NUMBER
        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        if( ransuu_tmp>0.5 ) then
           rho_up(i_count,i_count)=0.5+ransuu/2*((-1)**i_count)
        else
           rho_up(i_count,i_count)=0.5-ransuu/2*((-1)**i_count)
        end if

        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        if( ransuu_tmp<0.5 ) then
           rho_down(i_count,i_count)=0.5-ransuu/2*((-1)**(i_count-1))
        else
           rho_down(i_count,i_count)=0.5+ransuu/2*((-1)**(i_count-1))
        end if
     end do
  else if( HF_INITIAL_RHO==2 ) then            ! sdw的な解を仮定しない
     call fortran_random(ransuu,name,error)   ! 電子の存在確率は各サイトで
     call error_check(name,error)         ! 1になるようになっている。

     if( ransuu<0.5 ) then                ! ↑↓の順に並べる
        do i_count=1,TOTAL_SITE_NUMBER
           call fortran_random(ransuu,name,error)
           call error_check(name,error)

           rho_up(i_count,i_count)=ransuu
           rho_down(i_count,i_count)=1-ransuu
        end do
     else 
        do i_count=1,TOTAL_SITE_NUMBER    ! ↑↓の順
           call fortran_random(ransuu,name,error)
           call error_check(name,error)

           rho_up(i_count,i_count)=1-ransuu
           rho_down(i_count,i_count)=ransuu
        end do
     end if
  else if( HF_INITIAL_RHO==3 ) then       ! 全くランダムに作る(収束しにくい)
     do i_count=1,TOTAL_SITE_NUMBER
        call fortran_random(ransuu,name,error)
        call error_check(name,error) 
        rho_up(i_count,i_count)=ransuu
        
        call fortran_random(ransuu,name,error)
        call error_check(name,error)
        rho_down(i_count,i_count)=ransuu
     end do
  end if

! やっつけ
  if( JIGEN==2 .and. TOTAL_SITE_NUMBER==16 ) then
     rho_up(1,1)=1
     rho_up(2,2)=0
     rho_up(3,3)=1
     rho_up(4,4)=0

     rho_up(5,5)=0
     rho_up(6,6)=1
     rho_up(7,7)=0
     rho_up(8,8)=1

     rho_up(9,9)=1
     rho_up(10,10)=0
     rho_up(11,11)=1
     rho_up(12,12)=0

     rho_up(13,13)=0
     rho_up(14,14)=1
     rho_up(15,15)=0
     rho_up(16,16)=1
!
     rho_down(1,1)=0
     rho_down(2,2)=1
     rho_down(3,3)=0
     rho_down(4,4)=1

     rho_down(5,5)=1
     rho_down(6,6)=0
     rho_down(7,7)=1
     rho_down(8,8)=0

     rho_down(9,9)=0
     rho_down(10,10)=1
     rho_down(11,11)=0
     rho_down(12,12)=1

     rho_down(13,13)=1
     rho_down(14,14)=0
     rho_down(15,15)=1
     rho_down(16,16)=0
  end if

! 最大繰り返し回数まで計算する。それまでに収束したらループから抜ける。

  do i_count=1,MAX_HF_ITERATION
     call make_fock_matrix(rho_up,rho_down,name,error)
     call error_check(name,error)

! 対角化(実対称行列)
     lapack_work=0
     call dsyev(lapack_job,lapack_uplo,lapack_nnn,unitary_d_up,lapack_lda,&
          lapack_www_up,lapack_work,lapack_lwork,lapack_ifail)
     lapack_work=0
     call dsyev(lapack_job,lapack_uplo,lapack_nnn,unitary_d_down,lapack_lda,&
          lapack_www_down,lapack_work,lapack_lwork,lapack_ifail)

! 新しいrhoを求める
     call calc_new_rho(rho_up_new,rho_down_new,name,error)
     call error_check(name,error)

! 初期のrhoと新しいrhoとの比較(収束判定)
     
     result_up=0 ; result_down=0             ! 初期化(必須)

     do j_count=1,TOTAL_SITE_NUMBER
        result_up=result_up &
             +( rho_up_new(j_count,j_count)-rho_up(j_count,j_count) )**2
        result_down=result_down &
             +( rho_down_new(j_count,j_count)-rho_down(j_count,j_count) )**2
     end do

     result_up=(result_up/TOTAL_SITE_NUMBER)**0.5          ! 収束判定式
     result_down=(result_down/TOTAL_SITE_NUMBER)**0.5      !

! 以下の条件で収束として、hf解がもとまったこととする。
     if( result_up<dble(HF_CONVERGENCE_CONDITION) .and. &
          &result_down<dble(HF_CONVERGENCE_CONDITION) ) then 
        exit
     end if
     
! 収束してなかったら、新しく作ったrho_up_new,rho_down_newを元に、さらに
! 収束するまで計算を繰り返す。
     rho_up=rho_up_new ; rho_down=rho_down_new
     rho_up_new=0 ; rho_down_new=0
  end do

  if( i_count>MAX_HF_ITERATION ) then
     hf_iteration_result=i_count-1  
     write(*,*) "i_count,MAX_HF_ITERATION=",i_count,MAX_HF_ITERATION
     write(*,*) "rho_up=",rho_up
     write(*,*) "rho_down=",rho_down
     name="hf" ; error=-2 ; return            ! 収束しなかった場合エラーを出す
  else
     hf_iteration_result=i_count 
  end if

  deallocate(rho_up_new,stat=ier)
  call stat_check("rho_up_new","hf",2,ier)
  deallocate(rho_down_new,stat=ier)
  call stat_check("rho_down_new","hf",2,ier)

  deallocate(lapack_www_up,stat=ier)
  call stat_check("lapack_www_up","hf",2,ier)
  deallocate(lapack_www_down,stat=ier)
  call stat_check("lapack_www_down","hf",2,ier)
  deallocate(lapack_work,stat=ier)
  call stat_check("lapack_work","hf",2,ier)

!!! E_HF=<HF|H|HF>を求める

!
! トランスファー
!
  hf_energy_transfer=0
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        if( neighbor_table(i_count,j_count)==1 ) then
           hf_energy_transfer=hf_energy_transfer+dble(TRANSFER)&
                           *(rho_up(i_count,j_count)+rho_down(i_count,j_count))
        end if
     end do
  end do
!
! クーロン
!
  hf_energy_coulomb=0
  do i_count=1,TOTAL_SITE_NUMBER
     hf_energy_coulomb=hf_energy_coulomb+dble(COULOMB)&
                          *rho_up(i_count,i_count)*rho_down(i_count,i_count)
  end do


  hf_energy=hf_energy_transfer+hf_energy_coulomb

end subroutine hf



#ifdef DEBUG_HF
program main
  use global_variables
  implicit none

  character(32) :: name
  integer :: error
! MPI関連
  integer :: IERROR
! local
  integer,dimension(1) :: ic
  real(8) :: hf_energy
  integer :: hf_iteration_result
! init
  name="main" ; error=0
  unitary_d_up=0 ; unitary_d_down=0 
  hf_energy=0 ; hf_iteration_result=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(ic,name,error)
  call error_check(name,error)

  call hf(hf_energy,hf_iteration_result,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN,"TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON
  write(*,*) "TOTAL_DOWN_ELECTRON=",TOTAL_DOWN_ELECTRON
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*) "hf_energy=",hf_energy

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_HF */
