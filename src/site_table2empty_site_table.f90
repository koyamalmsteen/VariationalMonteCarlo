!
! $Id: site_table2empty_site_table.f90,v 1.11 2002/12/20 05:44:24 k-yukino Exp $
!
#include "parameter.h"

subroutine site_table2empty_site_table(total_sigma_electron,site_table,&
                                       empty_site_table,name,error)
  implicit none

  include "mpif.h"

  integer,intent(in) :: total_sigma_electron
  integer,dimension(TOTAL_SITE_NUMBER),intent(in) :: site_table
  integer,dimension(TOTAL_SITE_NUMBER-total_sigma_electron),&
                                      &intent(out) :: empty_site_table
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: i_count,j_count
! init
  empty_site_table=0 ; name="site_table2empty_site_table" ; error=0

! argument check
  if( total_sigma_electron<0 ) then
     error=-92 ; return
  else if( TOTAL_SITE_NUMBER<total_sigma_electron ) then
     error=-912 ; return
   end if
 
! set j_count=1

  j_count=1
  do i_count=1,TOTAL_SITE_NUMBER         
     if( site_table(i_count)==0 ) then   
        empty_site_table(j_count)=i_count
        j_count=j_count+1                
     end if
  end do

!
! 与えられたTOTAL_SITE_NUMBER-total_sigma_electron(つまりemptyサイトの数)
! と実際にチェックした電子配列要素の数とが等しいかのチェック
!
  if( TOTAL_SITE_NUMBER-total_sigma_electron/=j_count-1 ) then
     error=-1 ; return
  end if

end subroutine site_table2empty_site_table



#ifdef DEBUG_SITE_TABLE2EMPTY_SITE_TABLE
program main
  implicit none

  include "mpif.h"

  integer :: total_sigma_electron
  integer,dimension(TOTAL_SITE_NUMBER) :: site_table
  integer,dimension(:),allocatable :: empty_site_table
  character(32) :: name
  integer :: error
! local
  integer :: ier,i_count
! MPI関連
  integer :: IERROR
! init
  total_sigma_electron=0 ; site_table=0 ; name="main" ; error=0
! init [local]
  ier=0 ; i_count=0

  call MPI_INIT(IERROR)

  write(*,*) "total_sigma_electron >" ; read(*,*) total_sigma_electron
  
  allocate(empty_site_table(TOTAL_SITE_NUMBER-total_sigma_electron),stat=ier)
  call stat_check("empty_site_table","site_table2empty_site_table",1,ier)
  empty_site_table=0

  do i_count=1,TOTAL_SITE_NUMBER
     write(*,*) "i_count >" ; read(*,*) site_table(i_count)
  end do

  call site_table2empty_site_table(total_sigma_electron,site_table,&
                                   empty_site_table,name,error)
  call error_check(name,error)

  write(*,*) "JIGEN=",JIGEN
  write(*,*) "TOTAL_SITE_NUMBER=",TOTAL_SITE_NUMBER
  write(*,*) "empty_site_table=",empty_site_table

  call MPI_FINALIZE(IERROR)

end program main
#endif
