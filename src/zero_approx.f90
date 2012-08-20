!
! $Id: zero_approx.f90,v 1.6 2003/11/09 13:02:10 k-yukino Exp $
!
#include "parameter.h"

subroutine zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(out) :: zero_approx_energy
  integer,intent(out) :: hf_iteration_result 
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  real(8),dimension(:,:),allocatable :: tmp_unitary
  real(8) :: det_tmp_unitary
  integer :: i_count,j_count
  integer,dimension(:),allocatable :: ipiv
  integer :: ier,info
! init
  zero_approx_energy=0 ; hf_iteration_result=0
  name="zero_approx" ; error=0
! init local
  ier=0 ; info=0
  det_tmp_unitary=0

  if( INIT_WAVE_VECTOR==0 ) then
     call plane(zero_approx_energy,name,error)
     call error_check(name,error)
  else if( INIT_WAVE_VECTOR==1 ) then
     call hf(zero_approx_energy,hf_iteration_result,name,error)
     call error_check(name,error)
  else if ( INIT_WAVE_VECTOR==2 ) then
     call plane_with_field(zero_approx_energy,name,error)
     call error_check(name,error) 
  end if

!!
!! unitary$B9TNs$N(Bdeterminant$B$O(B1$B$K$J$k@-<A$rK\Ev$KK~$?$7$F$$$k$+3NG'$9$k!#(B
!!
  allocate(tmp_unitary(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("tmp_unitary","solve_hf_enenrgy",1,ier)
  tmp_unitary=0

  allocate(ipiv(TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("ipiv","main",1,ier)
  ipiv=0
! up
  tmp_unitary=unitary_d_up
! LU factorization
  ipiv=0
  call dgetrf(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER,&
              tmp_unitary,TOTAL_SITE_NUMBER,ipiv,info)
! determinant$B$r5a$a$k(B
  det_tmp_unitary=1
  do i_count=1,TOTAL_SITE_NUMBER
     det_tmp_unitary=det_tmp_unitary*tmp_unitary(i_count,i_count)
  end do

! $B$b$7(Bdet$B$,(B1$B$K6a$/$J$+$C$?$i%(%i!<(B
  if( abs(1-abs(det_tmp_unitary))>0.1 ) then 
     write(*,*) "det_tmp_unitary(up)=",det_tmp_unitary
     error=-1
     return
  end if

! down
  tmp_unitary=unitary_d_down

! LU$BJ,2r$7$F(B
  ipiv=0

  call dgetrf(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER,&
              tmp_unitary,TOTAL_SITE_NUMBER,ipiv,info)

! solve the determinant
  det_tmp_unitary=1
  do i_count=1,TOTAL_SITE_NUMBER
     det_tmp_unitary=det_tmp_unitary*tmp_unitary(i_count,i_count)
  end do
! $B$b$7(Bdet$B$,(B1$B$K6a$/$J$+$C$?$i%(%i!<(B
  if( abs(1-abs(det_tmp_unitary))>0.1 ) then
     write(*,*) "det_tmp_unitary(down)=",det_tmp_unitary
     error=-1
     return
  end if
  
  deallocate(ipiv,stat=ier)
  call stat_check("ipiv","main",2,ier)
  deallocate(tmp_unitary,stat=ier)
  call stat_check("tmp_unitary","main",2,ier)

!!
!! $B$5$F!"D69bB.2=$NE>CV$r$9$k!#(B($B$3$l0J9_%f%K%?%j!<9TNs$OE>CV$5$l$F$$$k(B)
!!
  allocate(tmp_unitary(TOTAL_SITE_NUMBER,TOTAL_SITE_NUMBER),stat=ier)
  call stat_check("tmp_unitary","solve_hf_enenrgy",1,ier)
! up
  tmp_unitary=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        tmp_unitary(j_count,i_count)=unitary_d_up(i_count,j_count)
     end do
  end do
  
  do i_count=1,TOTAL_SITE_NUMBER 
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=tmp_unitary(i_count,j_count)
     end do
  end do

! down
  tmp_unitary=0

  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        tmp_unitary(j_count,i_count)=unitary_d_down(i_count,j_count)
     end do
  end do
  
  do i_count=1,TOTAL_SITE_NUMBER 
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_down(i_count,j_count)=tmp_unitary(i_count,j_count)
     end do
  end do

  deallocate(tmp_unitary,stat=ier)
  call stat_check("tmp_unitary","solve_hf_energy",2,ier)

end subroutine zero_approx



#ifdef DEBUG_ZERO_APPROX
program main
  implicit none
  
  real(8) :: zero_approx_energy
  integer :: hf_iteration_result
  character(32) :: name
  integer :: error
! MPI
  integer :: IERROR
! local
  integer,dimension(1) :: ic
! init
  zero_approx_energy=0 ; hf_iteration_result=0 ; name="main" ; error=0
  ic(1)=1
  
  call MPI_INIT(IERROR)

! initialize parameters
  call init_global_variables(name,error)
  call error_check(name,error)

! initialize the seed of random number
  call random_init(ic,name,error)                 
  call error_check(name,error)

! zero approx
  call zero_approx(zero_approx_energy,hf_iteration_result,name,error)
  call error_check(name,error)
  
  write(*,*) "INIT_WAVE_VECTOR=",INIT_WAVE_VECTOR
  write(*,*) "JIGEN=",JIGEN,"TATAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON
  write(*,*) "TOTAL_DOWN_ELECTRON=",TOTAL_DOWN_ELECTRON
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*) "zero_approx_energy=",zero_approx_energy
  write(*,*) "hf_iteration_result=",hf_iteration_result

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_ZERO_APPROX */
