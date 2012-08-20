!
! $Id: invert.f90,v 1.6 2003/11/05 09:23:14 k-yukino Exp $
!
#include "parameter.h"

subroutine invert(spin,algorithm,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer,intent(in) :: spin,algorithm
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier
! BLAS
  integer,dimension(TOTAL_UP_ELECTRON) :: tmp_ipiv_up
  real(8),dimension(TOTAL_UP_ELECTRON) :: tmp_work_up
  integer,dimension(TOTAL_DOWN_ELECTRON) :: tmp_ipiv_down
  real(8),dimension(TOTAL_DOWN_ELECTRON) :: tmp_work_down
  integer :: info
! init
  name="invert" ; error=0 
  ier=0

!
!★１ LU分解とInvertで一連の作業である。↑スピン
!
  if( spin==1 ) then                           ! ↑スピン
     if( algorithm==1 ) then
        d_tilde_up_inverse=d_tilde_up
     else if( algorithm==2 ) then
        d_tilde_up_inverse=d_tilde_up_work
     else
!        d_tilde_up_inverse=d_tilde_up_for_ss
     end if

     ipiv_up=0

     tmp_ipiv_up=0
     tmp_work_up=0

! 以下のルーチンは、ノーマルlapackとnagと等価です。DEC alphやPentium III
! はATLAS利用で高速化が期待できるので、lapackを使用します。フリーだし。

! lapack(LU分解)
     call dgetrf(TOTAL_UP_ELECTRON,TOTAL_UP_ELECTRON,&
                 d_tilde_up_inverse,TOTAL_UP_ELECTRON,tmp_ipiv_up,info)

     after_lu_up=d_tilde_up_inverse
     ipiv_up=tmp_ipiv_up                 
     
! lapack(Invert)

     call dgetri(TOTAL_UP_ELECTRON,d_tilde_up_inverse,&
                 TOTAL_UP_ELECTRON,tmp_ipiv_up,tmp_work_up,&
                 TOTAL_UP_ELECTRON,info)
  else if( spin==2 ) then                                 ! ↓スピン
     if( algorithm==1 ) then 
        d_tilde_down_inverse=d_tilde_down
     else if( algorithm==2 ) then
        d_tilde_down_inverse=d_tilde_down_work
     else
!        d_tilde_down_inverse=d_tilde_down_for_ss
     end if

     ipiv_down=0

     tmp_ipiv_down=0
     tmp_work_down=0

! 以下のルーチンは、ノーマルlapackとnagと等価です。DEC alphやPentium III
! はATLAS利用で高速化が期待できるので、lapackを使用します。フリーだし。

! lapack(LU分解)
     call dgetrf(TOTAL_DOWN_ELECTRON,TOTAL_DOWN_ELECTRON,&
                 d_tilde_down_inverse,TOTAL_DOWN_ELECTRON,tmp_ipiv_down,info)
     after_lu_down=d_tilde_down_inverse
     ipiv_down=tmp_ipiv_down

! lapack(Invert)

     call dgetri(TOTAL_DOWN_ELECTRON,d_tilde_down_inverse,&
                 TOTAL_DOWN_ELECTRON,tmp_ipiv_down,tmp_work_down,&
                 TOTAL_DOWN_ELECTRON,info)
  end if

end subroutine invert



#ifdef DEBUG_INVERT
program main
  implicit none

  include "mpif.h"

  integer :: total_sigma_electron
  real(8),dimension(:,:),allocatable :: d_tilde_sigma,d_tilde_sigma_inverse
  character(32) :: name
  integer :: error
! MPI関連
  integer :: IERROR
! local
  integer :: ier,i_count,j_count
! init
  name="main" ; error=0

  call MPI_INIT(IERROR)

  write(*,*) "total_sigma_electron >" ; read(*,*) total_sigma_electron

  allocate(d_tilde_sigma(total_sigma_electron,total_sigma_electron),stat=ier)
  call stat_check("d_tilde_sigma","main",1,ier)
  allocate(d_tilde_sigma_invert(total_sigma_electron,total_sigma_electron),stat=ier)
  call stat_check("d_tilde_sigma","main",1,ier)
  
  do i_count=1,total_sigma_electron
     do j_count=1,total_sigma_electron
        write(*,*) "d_tilde_sigma >" ; read(*,*) d_tilde_sigma(i_count,j_count)
     end do
  end do

  call invert(total_sigma_electron,d_tilde_sigma,d_tilde_sigma_inverse,&
              name,error)
  call error_check(name,error)

  write(*,*) "d_tilde_sigma=",d_tilde_sigma
  write(*,*) "d_tilde_sigma_inverse=",d_tilde_sigma_inverse

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_INVERT */
