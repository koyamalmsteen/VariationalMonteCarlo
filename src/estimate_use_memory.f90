!
! $Id: estimate_use_memory.f90,v 1.2 2003/11/18 08:05:41 k-yukino Exp $
!
#include "parameter.h"

subroutine estimate_use_memory(result,name,error)
  use global_variables
  implicit none

  include "mpif.h"

  integer :: result
  character(32),intent(out) :: name
  integer,intent(out) :: error
! local
  integer :: four_byte,sixteen_byte
! init
  result=0 
  name="estimate_use_memory" ; error=0

  call est_global_variables(four_byte,sixteen_byte)
  call est_cor(four_byte,sixteen_byte)

  result=four_byte*4+sixteen_byte*16

end subroutine estimate_use_memory

subroutine est_global_variables(four_byte,sixteen_byte)
  use global_variables
  implicit none

  integer,intent(inout) :: four_byte,sixteen_byte

  four_byte=four_byte+TOTAL_SITE_NUMBER**2+TOTAL_SITE_NUMBER**2+4*TOTAL_SITE_NUMBER+TOTAL_SITE_NUMBER**2+1+2*TOTAL_SITE_NUMBER+2*TOTAL_SITE_NUMBER+TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON
  sixteen_byte=sixteen_byte+2*TOTAL_SITE_NUMBER**2+TOTAL_SITE_NUMBER+2*TOTAL_UP_ELECTRON**2+2*TOTAL_DOWN_ELECTRON**2+TOTAL_UP_ELECTRON**2+TOTAL_DOWN_ELECTRON**2+TOTAL_UP_ELECTRON**2+TOTAL_DOWN_ELECTRON**2+TOTAL_SITE_NUMBER+TOTAL_SITE_NUMBER
  /* 最後にlocal_vector1とlocal_vector2の分を足しておく */
  four_byte=four_byte+TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON+3*(MAX_MONTECARLO_SAMPLE-1)
  sixteen_byte=sixteen_byte+5+5*(MAX_MONTECARLO_SAMPLE-1)

end subroutine est_global_variables



subroutine est_cor(four_byte,sixteen_byte)
  use global_variables
  implicit none

  integer,intent(inout) :: four_byte,sixteen_byte

  four_byte=four_byte+1+TOTAL_UP_ELECTRON+TOTAL_DOWN_ELECTRON+3+3+4+3*7+3*8+4
!  sixteen_byte=sixteen_byte+1+3+4*MAX_MONTECARLO_SAMPLE*max_number_xi+MAX_MONTECARLO_SAMPLE*TOTAL_SITE_NUMBER+4*MAX_MONTECARLO_SAMPLE*max_number_xi+MAX_MONTECARLO_SAMPLE*TOTAL_SITE_NUMBER
  sixteen_byte=sixteen_byte+1+3

end subroutine est_cor
