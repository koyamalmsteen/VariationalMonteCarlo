!
! $Id: display_h_pow.f90,v 1.1 2003/11/18 08:05:41 k-yukino Exp $
!
subroutine display_h_pow(total_h0,total_h1,total_h2,total_h3,total_h4,&
     total_h5,name,error)
  implicit none

  real(8),intent(in) :: total_h0,total_h1,total_h2,total_h3,total_h4,total_h5
  character(32),intent(out) :: name
  integer,intent(out) :: error
! 
  name="display_h_pow" ; error=0
  
  write(*,*) "total_h0=",total_h0
  write(*,*) "total_h1=",total_h1
  write(*,*) "total_h2=",total_h2
  write(*,*) "total_h3=",total_h3

  write(*,*) "total_h1/total_h0=",total_h1/total_h0
  write(*,*) "total_h2/total_h1=",total_h2/total_h1,&
       "total_h3/total_h2=",total_h3/total_h2
  write(*,*) "total_h4=",total_h4
  write(*,*) "total_h5=",total_h5
  write(*,*) "total_h4/total_h3=",total_h4/total_h3,&
       "total_h5/total_h4=",total_h5/total_h4

end subroutine display_h_pow
