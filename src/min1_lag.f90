!
! $Id: min1_lag.f90,v 1.4 2004/02/15 10:52:35 k-yukino Exp ju448 $
!

#include "parameter.h"

subroutine min1_lag(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: h0,h1,h2,h3
  real(8),intent(out) :: pow1_c0,pow1_c1
  real(8),intent(out) :: result_pow1
  character(32),intent(out) :: name
  integer,intent(out) :: error
! scalapack
  integer :: ifail
  real(8),dimension(2) :: wr,wi
  real(8),dimension(8) :: work
  real(8),dimension(2,2) :: vl,vr
  real(8),dimension(2,2) :: aaa,bbb,ccc
  integer :: info
  integer,dimension(2) :: ipiv
! init
  pow1_c0=0 ; pow1_c1=0 ; result_pow1=0
  name="min1_lag" ; error=0
! init
  aaa=0 ; bbb=0 ; ccc=0
  ifail=0

  wr=0 ; wi=0
  vr=0 ; vl=0

  work=0

  aaa(1,1)=h0 ; aaa(1,2)=h1 
  aaa(2,1)=h1 ; aaa(2,2)=h2

  bbb(1,1)=h1 ; bbb(1,2)=h2 
  bbb(2,1)=h2 ; bbb(2,2)=h3

! aaaの逆行列を求める
  call dgetrf(2,2,aaa,2,ipiv,info)                    ! LU分解
  call dgetri(2,aaa,2,ipiv,work,2,info)               ! Invert

! aaaの逆行列とbbbとの行列積をcccとする

  ccc=0 

  ccc(1,1)=aaa(1,1)*bbb(1,1)+aaa(1,2)*bbb(2,1)
  ccc(1,2)=aaa(1,1)*bbb(1,2)+aaa(1,2)*bbb(2,2)
  ccc(2,1)=aaa(2,1)*bbb(1,1)+aaa(2,2)*bbb(2,1)
  ccc(2,2)=aaa(2,1)*bbb(1,2)+aaa(2,2)*bbb(2,2)

  ! 右固有ベクトルは求める(高校で習うやつね) 
  call dgeev('n','v',2,ccc,2,wr,wi,vl,2,vr,2,work,8,ifail)

  /*
   * wr(1)に対応するのがvr(1,1),vr(2,1)
   * wr(2)に対応するのがvr(1,2),vr(2,2)
   */

  if( wr(1)<wr(2) ) then
     result_pow1=wr(1)
     pow1_c0=vr(1,1)
     pow1_c1=vr(2,1)
  else
     result_pow1=wr(2)
     pow1_c0=vr(1,2)
     pow1_c1=vr(2,2)
  end if

end subroutine min1_lag



#ifdef DEBUG_MIN1_LAG
program main
  implicit none

  real(8) :: h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1
  character(32) :: name
  integer :: error
! init
  pow1_c1=0 ; pow1_c1=0 ; result_pow1=0

  h0=1000.00000000000
  h1=-4291.83682973980
  h2=38064.0073731103
  h3=-222327.386597429

  call min1_lag(h0,h1,h2,h3,pow1_c0,pow1_c1,result_pow1,name,error)

  write(*,*) "result_pow1=",result_pow1
  write(*,*) "pow1_c0=",pow1_c0,"pow1_c1=",pow1_c1

end program main
#endif /* DEBUG_MIN1_LAG */
