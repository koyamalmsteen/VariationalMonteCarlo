!
! $Id: choice_new_gamma_part1.f90,v 1.8 2002/12/20 05:44:30 k-yukino Exp $
!

! 電子を1個、ランダムに動かす。

#include "parameter.h"

subroutine choice_new_gamma_part1(gamma_up,gamma_down,idou_moto_electron_spin,&
                                  idou_moto_electron,idou_saki_site,name,error)
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: gamma_down
  integer,intent(out) :: idou_moto_electron_spin
  integer,intent(out) :: idou_moto_electron,idou_saki_site
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier,i_count
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer,dimension(:),allocatable :: empty_site_table
  integer :: number_occupy,total_empty_site_number,ransuu_int
! init 
  idou_moto_electron_spin=0 ; idou_moto_electron=0
  name="choice_new_gamma" ; error=0 ; ier=0 
! init[local]
  i_count=0
  total_empty_site_number=0 ; number_occupy=0
  site_table=0

!
! make a random number (1 to TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON)
! and make a choice electron will move.
!

! ↑と↓の電子をごっちゃにして、そこから移動元電子を一つ選び出す
  call fortran_random2(1,TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON,&
                      idou_moto_electron,name,error)
  call error_check(name,error)

! 1〜TOTAL_UP_ELECTRONまでなら↑の電子が選ばれる。
! TOTAL_UP_ELECTORN+1〜TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRONまでなら↓電子。
  if( idou_moto_electron<=TOTAL_UP_ELECTRON .and. idou_moto_electron>=1 ) then
     idou_moto_electron_spin=1
     idou_moto_electron=idou_moto_electron          ! このままでいい
  else if( idou_moto_electron>=TOTAL_UP_ELECTRON+1 .and. & 
           idou_moto_electron<=TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON ) then
     idou_moto_electron_spin=2
     idou_moto_electron=idou_moto_electron-TOTAL_UP_ELECTRON
  else
     error=-1                        ! 予期せぬエラー
  end if

! この段階で移動元は決まった。移動先はこれから空サイトを調べて、そこから
! ランダムに選ぶ。(サイトテーブルに変換して空サイトを調べる)
  if(  idou_moto_electron_spin==1 ) then
     do i_count=1,TOTAL_UP_ELECTRON
        site_table(gamma_up(i_count))=1
     end do
  else
     do i_count=1,TOTAL_DOWN_ELECTRON
        site_table(gamma_down(i_count))=1
     end do
  end if

! 空サイト数を調べる(ここではまず数だけ調べる)
  do i_count=1,TOTAL_SITE_NUMBER
     if( site_table(i_count)==1 ) then
        number_occupy=number_occupy+1
     else if( site_table(i_count)==0 ) then
        total_empty_site_number=total_empty_site_number+1
     else
        name="choice_gamma_sigma" ; error=-1 ; return
     end if
  end do

  if( total_empty_site_number==0 ) then         ! もし空サイトが全く無ければ、
     idou_moto_electron_spin=-1 ; return        
  end if

! 調べた空サイト数だけのメモリを確保して、空サイトのテーブルを作る。
  allocate(empty_site_table(total_empty_site_number),stat=ier)
  call stat_check("empty_site_table","choice_gamma_sigma",1,ier)

! 空サイトのテーブルを実際に作る。
  if( idou_moto_electron_spin==1 ) then
     call site_table2empty_site_table(TOTAL_UP_ELECTRON,site_table,&
                                      empty_site_table,name,error)
     call error_check(name,error)
  else
     call site_table2empty_site_table(TOTAL_DOWN_ELECTRON,site_table,&
                                      empty_site_table,name,error)
     call error_check(name,error)
  end if

! 移動先を探す
  call fortran_random2(1,total_empty_site_number,&
                      ransuu_int,name,error)
  call error_check(name,error)

  if( idou_moto_electron_spin==1 ) then
     gamma_up(idou_moto_electron)=empty_site_table(ransuu_int)
     idou_saki_site=empty_site_table(ransuu_int)
  else 
     gamma_down(idou_moto_electron)=empty_site_table(ransuu_int)
     idou_saki_site=empty_site_table(ransuu_int)
  end if

  deallocate(empty_site_table,stat=ier)
  call stat_check("empty_site_table","choice_gamma_sigma",2,ier)

end subroutine choice_new_gamma_part1



#ifdef DEBUG_CHOICE_NEW_GAMMA_PART1
program main
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count
  integer :: idou_moto_electron_spin,idou_moto_electron,idou_saki_site
! MPI関連
  integer :: IERROR 
! init
  name="main" ; error=0
! init [local]
  ier=0 ; i_count=0
  idou_moto_electron_spin=0 ; idou_moto_electron=0 ; idou_saki_site=0
! init [MPI] 
  IERROR=0

  call MPI_INIT(IERROR)

  call random_init(name,error)
  call error_check(name,error)

  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "gamma_up",i_count,">" ; read(*,*) gamma_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "gamma_down",i_count,">" ; read(*,*) gamma_down(i_count)
  end do

  call choice_new_gamma_part1(gamma_up,gamma_down,idou_moto_electron_spin,&
                              idou_moto_electron,idou_saki_site,name,error)
  call error_check(name,error)


  write(*,*) "JIGEN=",JIGEN
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "gamma_up=",gamma_up
  write(*,*) "gamma_down=",gamma_down

  write(*,*) "idou_moto_electron_spin=",idou_moto_electron_spin
  write(*,*) "idou_moto_electron=",idou_moto_electron
  write(*,*) "idou_saki_site=",idou_saki_site

  call MPI_FINALIZE(IERROR)

end program main

#endif /* DEBUG_CHOICE_NEW_GAMMA_PART1 */
