!
! $Id: choice_gamma_sigma_projection.f90,v 1.1 2004/02/15 10:52:35 k-yukino Exp $
!
#include "parameter.h"

subroutine choice_gamma_sigma_projection(total_sigma_electron,gamma_sigma,&
     name,error)
  use global_variables
  implicit none

  include "mpif.h"
  
  integer,intent(in) :: total_sigma_electron
  integer,dimension(total_sigma_electron),intent(inout) :: gamma_sigma
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer,dimension(total_sigma_electron) :: dummy     ! ダミー未使用
  integer :: i_count,j_count,k_count,l_count,ransuu_int,ier,&
             total_empty_site_number
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer,dimension(:),allocatable :: empty_site_table
! init arguments
  name="choice_gamma_sigma" ; error=0
! local init
  i_count=0 ; j_count=0 ; k_count=0 ; ransuu_int=0 ; ier=0
  site_table=0 ; total_empty_site_number=TOTAL_SITE_NUMBER

! argument check
  if( total_sigma_electron<=0 ) then       ! 電子数が0の時も対応していない
     error=-92 ; return
  else if( TOTAL_SITE_NUMBER<total_sigma_electron ) then
     error=-912 ; return
  end if

! これから初期電子状態を選ぶので、一つも電子は埋められていない
! だからサイトテーブルもemptyでよい

! 空サイトを保持しておく領域確保

  do i_count=1,total_sigma_electron
! 空のサイトを調べて、その中から一つのサイトを電子で埋める
     allocate(empty_site_table(total_empty_site_number),stat=ier)
     call stat_check("empty_site_table","choice_gamma_sigma",1,ier)
     empty_site_table=0  

     total_empty_site_number=0
     do l_count=1,TOTAL_SITE_NUMBER
        if( site_table(l_count)==0 ) then
           total_empty_site_number=total_empty_site_number+1
        end if
     end do

!  以下のようにして空サイトのチェックを行なう。

     k_count=1

     do j_count=1,TOTAL_SITE_NUMBER
        if( site_table(j_count)==0 ) then
           empty_site_table(k_count)=j_count
           k_count=k_count+1
        end if
     end do

! チェックした空サイトから電子を埋めるサイトをひとつ決める

     call fortran_random2(1,total_empty_site_number,ransuu_int,name,error)
     call error_check(name,error)

     gamma_sigma(i_count)=empty_site_table(ransuu_int)

!
! gamma_sigmaがひとつ電子で埋められたので、新しいsite_tableを作る必要ができた。
! つまり、あたらしいsite_tableを作り、上記の空サイトを調べて、新しく電子を
! 詰める作業をするわけである。
!  ここでは、以下にサイトへの変換プログラムを書く。
!
     deallocate(empty_site_table,stat=ier)
     call stat_check("empty_site_table","choice_gamma_sigma",2,ier)

! サイト表示への変換
     site_table=0

     do j_count=1,i_count
        site_table(gamma_sigma(j_count))=1
     end do

     total_empty_site_number=total_empty_site_number-1
  end do

! slatec library(昇べき(1)に並べる。汎用性のため、cxmlからslatecに書き換えた)
  call isort(gamma_sigma,dummy,total_sigma_electron,1)

end subroutine choice_gamma_sigma_projection



#ifdef DEBUG_CHOICE_GAMMA_SIGMA
program main
  implicit none

  include "mpif.h"

  integer :: total_sigma_electron
  integer,dimension(:),allocatable :: gamma_sigma
  character(32) :: name
  integer :: error
! local
  integer :: i_count,ier
  integer,external :: iargc
  character(32) :: string
! MPI関連
  integer :: IERROR
! init
  name="main" ; error=0
  i_count=0 ; ier=0

  call MPI_INIT(IERROR)

  if( iargc()==2 )then
     call getarg(2,string) ; read(string,*) total_sigma_electron
  else
     write(*,*) "total_sigma_electron >" ; read(*,*) total_sigma_electron
  end if

  allocate(gamma_sigma(total_sigma_electron),stat=ier)
  call stat_check("gamma_sigma","choice_gamma_sigma",1,ier)

! init
  gamma_sigma=0

! 乱数テーブルの初期化 (要注意箇所)
  call random_init(name,error)
  call error_check(name,error)

  call choice_gamma_sigma(total_sigma_electron,gamma_sigma,name,error)
  call error_check(name,error)

  write(*,*) "gamma_sigma=",gamma_sigma

!
  deallocate(gamma_sigma,stat=ier)
  call stat_check("gamma_sigma","choice_gamma_sigma",2,ier)

  
  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_CHOICE_GAMMA_SIGMA */ 
