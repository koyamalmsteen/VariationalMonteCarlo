!
! $Id$
!
#include "parameter.h"

subroutine plane_with_field(plane_energy,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  real(8),intent(out) :: plane_energy
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
  integer :: ier
  real(8) :: plane_energy_transfer,plane_energy_coulomb
! scalapack
  integer :: nnn,lda,lwork,ifail
  real(8),dimension(:),allocatable :: www_up,www_down
  real(8),dimension(:),allocatable :: work
  character*1 :: job,uplo
! init
  plane_energy=0
  name="plane" ; error=0
! init local
  ier=0 ; plane_energy_transfer=0 ; plane_energy_coulomb=0

  job='v'                   ! eigen value and eigen vector are needed
  uplo='u'                  ! upper triangular part of A is stored
  nnn=TOTAL_SITE_NUMBER
  lda=TOTAL_SITE_NUMBER
  lwork=64*TOTAL_SITE_NUMBER
  ifail=0

  allocate(www_up(lda),stat=ier)
  call stat_check("www_up","plane",1,ier)
  allocate(www_down(lda),stat=ier)
  call stat_check("www_down","plane",1,ier)
  allocate(work(lwork),stat=ier)
  call stat_check("work","plane",1,ier)
  www_up=0 ; www_down=0 ; work=0

! make fock matrix
  do i_count=1,TOTAL_SITE_NUMBER
     do j_count=1,TOTAL_SITE_NUMBER
        unitary_d_up(i_count,j_count)=dble(TRANSFER) &
             *neighbor_table(i_count,j_count)
        unitary_d_down(i_count,j_count)=dble(TRANSFER) &
             *neighbor_table(i_count,j_count)

        if( i_count==j_count ) then
           unitary_d_up(i_count,j_count)=unitary_d_up(i_count,j_count) &
                -electric_field*xaxis(i_count)
           unitary_d_down(i_count,j_count)=unitary_d_down(i_count,j_count) &
                -electric_field*xaxis(i_count)
        end if
     end do
  end do

! 対角化(実対称行列)
  work=0
  call dsyev(job,uplo,nnn,unitary_d_up,lda,www_up,work,lwork,ifail)
  work=0
  call dsyev(job,uplo,nnn,unitary_d_down,lda,www_down,work,lwork,ifail)

  plane_energy_transfer=0

  do i_count=1,TOTAL_UP_ELECTRON
     plane_energy_transfer=plane_energy_transfer+www_up(i_count)
  end do

  do i_count=1,TOTAL_DOWN_ELECTRON
     plane_energy_transfer=plane_energy_transfer+www_down(i_count)
  end do
  
  plane_energy_coulomb=0

  plane_energy=plane_energy_transfer+dble(TOTAL_SITE_NUMBER)*dble(COULOMB)/dble(4.0)

  deallocate(www_up,stat=ier)
  call stat_check("www_up","plane",2,ier)
  deallocate(www_down,stat=ier)
  call stat_check("www_down","plane",2,ier)
  deallocate(work,stat=ier)
  call stat_check("work","plane",2,ier)

end subroutine plane_with_field



#ifdef DEBUG_PLANE
program main
  use global_variables
  implicit none

  character(32) :: name
  integer :: error
! local
  integer,dimension(2) :: ic
  real(8) :: plane_energy
! 
  name="main" ; error=0 
  ic=1 ; plane_energy=0

  call MPI_INIT(IERROR)

! パラメーターの初期化
  call init_global_variables(name,error)
  call error_check(name,error)

! 乱数テーブルの初期化
  call random_init(ic,name,error)
  call error_check(name,error)

! 平面波解
  call plane(plane_energy,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN,"TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "TOTAL_UP_ELECTRON=",TOTAL_UP_ELECTRON
  write(*,*) "TOTAL_DOWN_ELECTRON=",TOTAL_DOWN_ELECTRON
  write(*,*) "COULOMB=",COULOMB,"TRANSFER=",TRANSFER
  write(*,*) "plane_energy=",plane_energy

  call MPI_FINALIZE(IERROR)

end program main
#endif /* DEBUG_PLANE */
