!
! $Id: allocate_global_variables.f90,v 1.1 2003/11/18 08:05:41 k-yukino Exp $
!
#include "parameter.h"

subroutine allocate_global_variables(name,error)
  use global_variables
  implicit none

  include "mpif.h"

  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: ier
! init
  name="allocate_global_variables" ; error=0
 
!
! $BAj4X4X?t(B
!

!
! sz
!
  allocate(result_sz0(max_number_xi),stat=ier)
  call stat_check("result_sz0","allocate_global_variables",1,ier)
  result_sz0=0

  allocate(result_sz1(max_number_xi),stat=ier)
  call stat_check("result_sz1","allocate_global_variables",name,error)
  result_sz1=0

  allocate(result_sz2(max_number_xi),stat=ier)
  call stat_check("result_sz2","allocate_global_variables",name,error)
  result_sz2=0
!
! sxsy
!
  allocate(result_sxsy0(max_number_xi),stat=ier)
  call stat_check("result_sxsy0","allocate_global_variables",1,ier)
  result_sxsy0=0

  allocate(result_sxsy1(max_number_xi),stat=ier)
  call stat_check("result_sxsy1","allocate_global_variables",1,error)
  result_sxsy1=0

  allocate(result_sxsy2(max_number_xi),stat=ier)
  call stat_check("result_sxsy2","allocate_global_variables",1,error)
  result_sxsy2=0
!
! s
!
  allocate(result_s0(max_number_xi),stat=ier)
  call stat_check("result_s0","allocate_global_variables",1,ier)
  result_s0=0

  allocate(result_s1(max_number_xi),stat=ier)
  call stat_check("result_s1","allocate_global_variables",1,ier)
  result_s1=0

  allocate(result_s2(max_number_xi),stat=ier)
  call stat_check("result_s2","allocate_global_variables",1,ier)
  result_s2=0
!
! charge
!
  allocate(result_charge0(max_number_xi),stat=ier)
  call stat_check("result_charge0","allocate_global_variables",1,ier)
  result_charge0=0

  allocate(result_charge1(max_number_xi),stat=ier) 
  call stat_check("result_charge1","allocate_global_variables",name,error)
  result_charge1=0

  allocate(result_charge2(max_number_xi),stat=ier) 
  call stat_check("result_charge2","allocate_global_variables",name,error)
  result_charge2=0
!
! super
!
  allocate(result_super0(max_number_xi),stat=ier)
  call stat_check("result_super0","allocate_global_variables",1,ier)
  result_super0=0

  allocate(result_super1(max_number_xi),stat=ier)
  call stat_check("result_super1","allocate_global_variables",1,error)
  result_super1=0

  allocate(result_super2(max_number_xi),stat=ier) 
  call stat_check("result_super2","allocate_global_variables",1,error)
  result_super2=0

! result_charge_density0,1,2$B$O%5%$%H?t$HF1$8Bg$-$5$r;}$D$3$H$,J,$+$C$F$$$k$N$G(B
! $B@EE*$K3NJ]$7$F$$$k!#$G$b=i4|2=$O$3$3$G$7$F$"$2$k!#(B
  result_charge_density0=0
  result_charge_density1=0
  result_charge_density2=0

end subroutine allocate_global_variables
