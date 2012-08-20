! 
! $Id: choice_new_gamma_part2.f90,v 1.11 2002/12/20 06:15:13 k-yukino Exp k-yukino $
!

! 再近接サイトに、電子を動かす。再近接が全て埋まっていたら、
! idou_moto_electron_spin=-1を返す。呼び出し元で、それをチェックして
! 再度このルーチンを呼び出す。

#include "parameter.h"

subroutine choice_new_gamma_part2(gamma_up,gamma_down,idou_moto_electron_spin,&
                                  idou_moto_electron,idou_saki_site,name,error)
  use global_variables
  implicit none

  include "mpif.h"

! ちなみにidou_moto_electronはサイト表示ではないよ。ガンマ表示ね♪
! idou_saki_siteはサイト表示
  integer,intent(out) :: idou_moto_electron_spin,idou_moto_electron ,&
                         idou_saki_site 
  integer,dimension(TOTAL_UP_ELECTRON),intent(inout) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON),intent(inout) :: gamma_down
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier,ransuu_int
  integer,dimension(JIGEN*2) :: nearest_neighbor
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer :: i_count,j_count,temp
  integer :: empty_site_number
  integer,external :: state_include
  integer,dimension(:),allocatable :: empty_site_table
! init
  idou_moto_electron_spin=0 ; idou_moto_electron=0
  idou_saki_site=0
  name="choice_new_gamma_part2" ; error=0
! local init
  i_count=0 ; site_table=0 ; temp=0

  do while (1)

! ↑と↓の電子をごっちゃにして、そこから移動元電子を一つ選び出す
     call fortran_random2(1,TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON,&
                          idou_moto_electron,name,error)
     call error_check(name,error)

! 1〜TOTAL_UP_ELECTRONまでなら↑の電子が選ばれる。
! TOTAL_UP_ELECTORN+1〜TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRONまでなら↓電子。

     if( idou_moto_electron>=1 .and. &
         idou_moto_electron<=TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON ) then

        if( idou_moto_electron<=TOTAL_UP_ELECTRON ) then

           idou_moto_electron_spin=1               ! upスピン 
           idou_moto_electron=idou_moto_electron   ! 1〜TOTAL_UP_ELECTRONまで
                                                   ! なのでそのまま 
        else if( idou_moto_electron>=TOTAL_UP_ELECTRON+1 ) then
           idou_moto_electron_spin=2               ! downスピン
           idou_moto_electron=idou_moto_electron  &! ここで、左記のように
                              -TOTAL_UP_ELECTRON   !  置き換える
        end if
     end if

!
! この段階で移動元は決まった。移動先はこの移動元に隣接する(JIGEN*2)サイトから
! ランダムに選ぶ。再近接サイトテーブルを作る。
!

     j_count=1                                   ! カウンタ  

     if( idou_moto_electron_spin==1 ) then       ! upスピンが選択された処理
        do i_count=1,TOTAL_SITE_NUMBER
           if( neighbor_table(gamma_up(idou_moto_electron),i_count)==1 ) then
              nearest_neighbor(j_count)=i_count
              j_count=j_count+1 
           end if
        end do
     else                                        ! downスピンが選択された処理
        do i_count=1,TOTAL_SITE_NUMBER
           if( neighbor_table(gamma_down(idou_moto_electron),i_count)==1 ) then
              nearest_neighbor(j_count)=i_count
              j_count=j_count+1
           end if
        end do
     end if

!
! 再近接サイトを網羅したテーブルができた。その中から、空サイトのみのテーブル
! を作り出す。
!
     empty_site_number=0

     if( idou_moto_electron_spin==1 ) then         ! upスピンが選択された処理
! まず、再近接空サイト数を調べる
        do i_count=1,JIGEN*2
! moshi hashi de nakereba
           if( nearest_neighbor(i_count)/=0 ) then
              if( state_include(TOTAL_UP_ELECTRON,&
                   &nearest_neighbor(i_count),gamma_up)/=1 ) then
                 empty_site_number=empty_site_number+1
              end if
           end if
        end do
     else                                          ! downスピンが選択された処理
! まず、再近接空サイト数を調べる
        do i_count=1,JIGEN*2
! moshi hashi de nakereba
           if( nearest_neighbor(i_count)/=0 ) then
              if( state_include(TOTAL_DOWN_ELECTRON,&
                   &nearest_neighbor(i_count),gamma_down)/=1 ) then
                 empty_site_number=empty_site_number+1
              end if
           end if
        end do
     end if

     if( empty_site_number==0 ) then            ! もし再近接が全て埋まって
        idou_moto_electron_spin=-1 ; return     ! いたら
     end if

! もし、移動元の再近接が全て、埋められていたら、以下の処理を行なわないで、
! 移動元を選び直し、始めからやりなおす。

     if( empty_site_number/=0 ) then               ! 移動先があれば
  
        if( idou_moto_electron_spin==1 ) then      ! upスピンが選択された処理
           allocate(empty_site_table(empty_site_number),stat=ier)
           call stat_check("empty_site_table","choice_new_gamma_part2",1,ier)
        
! 実際に再近接空サイトのテーブルを作る
           
           j_count=1
        
           do i_count=1,JIGEN*2
! moshi hashi denai nara
              if( nearest_neighbor(i_count)/=0 ) then
                 if( state_include(TOTAL_UP_ELECTRON,&
                      &nearest_neighbor(i_count),gamma_up)/=1 ) then
                    empty_site_table(j_count)=nearest_neighbor(i_count)
                    j_count=j_count+1
                 end if
              end if
           end do
        else                                    ! downスピンが選択された処理
           allocate(empty_site_table(empty_site_number),stat=ier)
           call stat_check("empty_site_table","choice_new_gamma_part2",1,ier)
        
! 実際に再近接空サイトのテーブルを作る

           j_count=1
        
           do i_count=1,JIGEN*2
! moshi hashi de nainara
              if( nearest_neighbor(i_count)/=0 ) then
                 if( state_include(TOTAL_DOWN_ELECTRON,&
                      &nearest_neighbor(i_count),gamma_down)/=1 ) then
                    empty_site_table(j_count)=nearest_neighbor(i_count)
                    j_count=j_count+1
                 end if
              end if
           end do
        end if

!以下の乱数で飛ばす先を決める。
        call fortran_random2(1,empty_site_number,ransuu_int,name,error)
        call error_check(name,error)

! では実際に置き換える(ここはサイト表示の方が便利なので、表示を変更したい
! 気持ちは山々だが、実は不都合があるのでしない)←アクセプタンスレシオを
! 求める高速化ルーチンの時困ります。
        if( idou_moto_electron_spin==1 ) then      ! upスピンが選択された処理
           gamma_up(idou_moto_electron)=empty_site_table(ransuu_int)
           idou_saki_site=empty_site_table(ransuu_int)
           return
        else                                       ! downスピンが選択された処理
           gamma_down(idou_moto_electron)=empty_site_table(ransuu_int)
           idou_saki_site=empty_site_table(ransuu_int)
           return
        end if
     end if
  end do

end subroutine choice_new_gamma_part2



#ifdef DEBUG_CHOICE_NEW_GAMMA_PART2
program main
  use global_variables
  implicit none

  include "mpif.h"

  integer,dimension(TOTAL_UP_ELECTRON) :: gamma_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: gamma_down
  character(32) :: name
  integer :: error
! local
  integer :: idou_moto_electron_spin,idou_moto_electron,idou_saki_site
  integer :: ier,i_count
! MPI関連
  integer :: IERROR
! init
  gamma_up=0 ; gamma_down=0 ; name="main" ; error=0
! init [local]
  idou_moto_electron_spin=0 ; idou_moto_electron=0 ; idou_saki_site=0
  ier=0 ; i_count=0
! init [MPI]
  IERROR=0

  call MPI_INIT(IERROR)

  call init_global_variables(name,error)
  call error_check(name,error)

  call random_init(name,error)
  call error_check(name,error)
 
  do i_count=1,TOTAL_UP_ELECTRON
     write(*,*) "gamma_up",i_count,">" ; read(*,*) gamma_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     write(*,*) "gamma_down",i_count,">" ; read(*,*) gamma_down(i_count)
  end do

  call choice_new_gamma_part2(gamma_up,gamma_down,&
                              idou_moto_electron_spin,idou_moto_electron,&
                              idou_saki_site,name,error)
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

#endif /* DEBUG_CHOICE_NEW_GAMMA_PART2 */
