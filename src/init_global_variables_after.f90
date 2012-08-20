!
! $Id: init_global_variables_after.f90,v 1.1 2004/03/10 03:04:02 k-yukino Exp $
!

#include "parameter.h"

subroutine init_global_variables_after(name,error)
  use global_variables
  implicit none

  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier
! init
  name="init_global_variables_after" ; error=0

!
! エラーバー(スピンSz)
!
  allocate(buff_cor5(5,max_number_xi),stat=ier) ; buff_cor5=0
! sz
  allocate(splited_sz0(5,max_number_xi),stat=ier) ; splited_sz0=0
  allocate(splited_sz1(5,max_number_xi),stat=ier) ; splited_sz1=0
  allocate(splited_sz2(5,max_number_xi),stat=ier) ; splited_sz2=0
! charge
  allocate(splited_charge0(5,max_number_xi),stat=ier) ; splited_charge0=0
  allocate(splited_charge1(5,max_number_xi),stat=ier) ; splited_charge1=0
  allocate(splited_charge2(5,max_number_xi),stat=ier) ; splited_charge2=0
! charge_density
  allocate(splited_charge_density0(5,TOTAL_SITE_NUMBER),stat=ier)
  splited_charge_density0=0
  allocate(splited_charge_density1(5,TOTAL_SITE_NUMBER),stat=ier)
  splited_charge_density1=0
  allocate(splited_charge_density2(5,TOTAL_SITE_NUMBER),stat=ier)
  splited_charge_density2=0
! sxsy
  allocate(splited_sxsy0(5,max_number_xi),stat=ier) ; splited_sxsy0=0
  allocate(splited_sxsy1(5,max_number_xi),stat=ier) ; splited_sxsy1=0
  allocate(splited_sxsy2(5,max_number_xi),stat=ier) ; splited_sxsy2=0
! s
  allocate(splited_s0(5,max_number_xi),stat=ier) ; splited_s0=0
  allocate(splited_s1(5,max_number_xi),stat=ier) ; splited_s1=0
  allocate(splited_s2(5,max_number_xi),stat=ier) ; splited_s2=0
! super
  allocate(splited_super0(5,max_number_xi),stat=ier) ; splited_super0=0
  allocate(splited_super1(5,max_number_xi),stat=ier) ; splited_super1=0
  allocate(splited_super2(5,max_number_xi),stat=ier) ; splited_super2=0
! pole(静的に確保済み)
  splited_pole0=0 ; splited_pole1=0 ; splited_pole2=0

end subroutine init_global_variables_after
