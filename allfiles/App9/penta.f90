!************************** P E N T A ***************************
!                                                               *
!     Solution of a linear system of algebraic equations with   *
!     a pentadiagonal matrix of coefficients.(No pivoting)      *
!     Equation no. i :                                          *
!                                                               *
!     e(i)*x(i-2)+ a(i)*x(i-1) + b(i)*x(i) +                    *
!     c(i)*x(i+1) + f(i)*x(i+2)  = d(i), i = 1,2,...n           *
!                                                               *
!                       === Use ===                             *
!                                                               *
!                  call penta(n,a,b,c,d,x)                      *
!                            or                                 *
!                  call penta(n,a,b,c,d,d)                      *
!                                                               *
!      In the last case, vector d contains the solution.        *
!                                                               *
!                     === Input ===                             *
!                                                               *
!     n ....... integer     . Number of equations               *
!     e(1:n) .. real vector . Lowest diagonal.Elements e(1)     *
!                             and e(2) are not used             *
!     a(1:n) .. real vector . Lower diagonal.Element a(1)       *
!                                            is not used.       *
!     b(1:n) .. real vector . Main diagonal                     *
!     c(1:n) .. real vector . Upper diagonal.Element c(n)       *
!                                            is not used.       *
!     f(1:n) .. real vector . Uppermost diagonal.Elements       *
!                             f(n-1)and f(n) are not used.      *
!     d(1:n) .. real vector . Right hand side of the system.    *
!                                                               *
!                     === Output ===                            *
!                                                               *
!     x(1:n) .. real vector . The solution vector               *
!                                                               *
!********************** fortran 90 ******************************

      subroutine penta(n,e,a,b,c,f,d,x)
	  implicit none
      integer, intent(in) :: n
      real, intent(inout), dimension(n) :: e,a,b,c,f,d
	  real, intent(out) :: x(n)
	  !  --- Local variables ---
	  integer :: i
	  real :: q
      do i = 2,n-1
         q = a(i)/b(i-1)
         b(i) = b(i) - q*c(i-1)
         c(i) = c(i) - q*f(i-1)
         d(i) = d(i) - q*d(i-1)
         q = e(i+1)/b(i-1)
         a(i+1) = a(i+1) - q*c(i-1)
         b(i+1) = b(i+1) - q*f(i-1)
         d(i+1) = d(i+1) - q*d(i-1)
      end do
      q = a(n)/b(n-1)
      b(n) = b(n) - q*c(n-1)
      x(n) = (d(n) - q*d(n-1))/b(n)
      x(n-1) = (d(n-1) - c(n-1)*x(n))/b(n-1)
      do i = n-2,1,-1
         x(i) = (d(i) - f(i)*x(i+2) - c(i)*x(i+1))/b(i)
      end do
      return
      end
