!
! $Id: min2_lag.f90,v 1.3 2004/02/15 10:52:35 k-yukino Exp ju448 $
!

#include "parameter.h"

subroutine min2_lag(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,pow2_c2,&
     result_pow2,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(in) :: h0,h1,h2,h3,h4,h5
  real(8),intent(out) :: pow2_c0,pow2_c1,pow2_c2
  real(8),intent(out) :: result_pow2
  character(32),intent(out) :: name
  integer,intent(out) :: error
!
  integer :: ifail
  real(8),dimension(3) :: wr,wi
  real(8),dimension(12) :: work
  character*1 :: jobvl,jobvr

  real(8),dimension(3,3) :: vl,vr
  real(8),dimension(3,3) :: aaa,bbb,ccc
  integer :: info
  integer,dimension(3) :: ipiv
! 
  integer :: i_count,j_count
!
  pow2_c0=0 ; pow2_c1=0 ; pow2_c2=0 
  result_pow2=0
  name="min2_lag" ; error=0

  aaa=0 ; bbb=0 ; ccc=0

  ifail=0

  wr=0 ; wi=0

  vr=0 ; vl=0

  work=0

  aaa(1,1)=h0 ; aaa(1,2)=h1 ; aaa(1,3)=h2
  aaa(2,1)=h1 ; aaa(2,2)=h2 ; aaa(2,3)=h3
  aaa(3,1)=h2 ; aaa(3,2)=h3 ; aaa(3,3)=h4

  bbb(1,1)=h1 ; bbb(1,2)=h2 ; bbb(1,3)=h3
  bbb(2,1)=h2 ; bbb(2,2)=h3 ; bbb(2,3)=h4
  bbb(3,1)=h3 ; bbb(3,2)=h4 ; bbb(3,3)=h5

! aaaの逆行列を求める(LU分解->invert)

  call dgetrf(3,3,aaa,3,ipiv,info)                ! LU分解
  call dgetri(3,aaa,3,ipiv,work,3,info)           ! Invert

! aaaの逆行列とbbbとの行列積をcccとする

  ccc=0

  ccc(1,1)=aaa(1,1)*bbb(1,1)+aaa(1,2)*bbb(2,1)+aaa(1,3)*bbb(3,1)
  ccc(1,2)=aaa(1,1)*bbb(1,2)+aaa(1,2)*bbb(2,2)+aaa(1,3)*bbb(3,2)
  ccc(1,3)=aaa(1,1)*bbb(1,3)+aaa(1,2)*bbb(2,3)+aaa(1,3)*bbb(3,3)
  ccc(2,1)=aaa(2,1)*bbb(1,1)+aaa(2,2)*bbb(2,1)+aaa(2,3)*bbb(3,1)
  ccc(2,2)=aaa(2,1)*bbb(1,2)+aaa(2,2)*bbb(2,2)+aaa(2,3)*bbb(3,2)
  ccc(2,3)=aaa(2,1)*bbb(1,3)+aaa(2,2)*bbb(2,3)+aaa(2,3)*bbb(3,3)
  ccc(3,1)=aaa(3,1)*bbb(1,1)+aaa(3,2)*bbb(2,1)+aaa(3,3)*bbb(3,1)
  ccc(3,2)=aaa(3,1)*bbb(1,2)+aaa(3,2)*bbb(2,2)+aaa(3,3)*bbb(3,2)
  ccc(3,3)=aaa(3,1)*bbb(1,3)+aaa(3,2)*bbb(2,3)+aaa(3,3)*bbb(3,3)

  ! 右固有ベクトルのみを求める 
  call dgeev('n','v',3,ccc,3,wr,wi,vl,3,vr,3,work,12,ifail)

  /*
   * wr(1)に対応するのがvr(1,1),vr(2,1),vr(3,1)
   * wr(2)に対応するのがvr(1,2),vr(2,2),vr(3,2)
   * wr(3)に対応するのがvr(1,3),vr(2,3),vr(3,3)
   */


  if( wr(1)<wr(2) ) then
     if( wr(1)<wr(3) ) then
        result_pow2=wr(1)
        pow2_c0=vr(1,1)
        pow2_c1=vr(2,1)
        pow2_c2=vr(3,1)
     else
        result_pow2=wr(3)
        pow2_c0=vr(1,3)
        pow2_c1=vr(2,3)
        pow2_c2=vr(3,3)
     end if
  else
     if( wr(2)<wr(3) ) then
        result_pow2=wr(2)
        pow2_c0=vr(1,2)
        pow2_c1=vr(2,2)
        pow2_c2=vr(3,2)
     else
        result_pow2=wr(3)
        pow2_c0=vr(1,3)
        pow2_c1=vr(2,3)
        pow2_c2=vr(3,3)
     end if
  end if

end subroutine min2_lag



#ifdef DEBUG_MIN2_LAG 
program main
  implicit none

  real(8) :: h0,h1,h2,h3,h4,h5
  real(8) :: pow2_c0,pow2_c1,pow2_c2,result_pow2
  character(32) :: name
  integer :: error
! init
  name="main" ; error=0 
  pow2_c0=0 ; pow2_c1=0 ; pow2_c2=0 ; result_pow2=0

  h0=1000.00000000000
  h1=-4291.83682973980
  h2=38064.0073731103
  h3=-222327.386597429
  h4=2263343.83085554
  h5=-10946209.6183627

  call  min2_lag(h0,h1,h2,h3,h4,h5,pow2_c0,pow2_c1,pow2_c2,&
       result_pow2,name,error)

  write(*,*) "result_pow2=",result_pow2
  write(*,*) "pow2_c0=",pow2_c0,"pow2_c1=",pow2_c1,"pow2_c2=",pow2_c2

end program main
#endif /* DEBUG_MIN2_LAG */
