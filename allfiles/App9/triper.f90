!**************************** T R I P E R ***************************
!                                                                   *
!     Solution of a tridiagonal system of equations with            *
!     periodicity in the first point ,point no.1, and last          *
!     point , point no. n + 1, giving n equations for the unknowns. *
!                                                                   *
!     a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)  ,k= 1,2,,,n     *
!                                                                   *
!     The element in the upper right corner is stored in a(1)       *
!     The element in the lower left  corner is stored in c(N)       *
!     a,b,c,d,x are vectors with dimension n ,to be                 *
!     specified in the calling program.                             *
!                                                                   *
!     Reference: C. Hirsch ," Numerical Computation of Internal     *
!                             and External Flows", Volume 1         *
!                             John Wiley & Sons  1988               *
!                                                                   *
!*************************** for90 **********************************

      subroutine triper(n,a,b,c,d,x)
	  implicit none
      integer, intent(in) :: n
      real, intent(inout), dimension(n) :: a,b,c,d
	  real, intent(out) :: x(n)
	  !  --- Local variables ---
	  integer :: k
	  real :: xn

      b(1) = 1.0/b(1)
      x(1) = -a(1)*b(1)
      a(1) = d(1)*b(1)
      ! --- Elimination ---
      do k = 2,n-1
         c(k-1) = c(k-1)*b(k-1)
         b(k) = b(k) - a(k)*c(k-1)
         b(k) = 1.0/b(k)
         x(k) = - a(k)*x(k-1)*b(k)
         a(k) = (d(k) - a(k)*a(k-1))*b(k)
      end do
      x(n-1) = x(n-1) - c(n-1)*b(n-1)
      ! --- Backsubstitution ---
      d(n-1) = a(n-1)
      b(n-1) = x(n-1)
      do k = n-2,1,-1
         d(k) = a(k) - c(k)*d(k+1)
         b(k) = x(k) - c(k)*b(k+1)
      end do
      xn = d(n) - c(n)*d(1) - a(n)*d(n-1)
      xn = xn/(b(n) + a(n)*b(n-1) + c(n)*b(1))
      x(n) = xn
      do k = 1,n-1
         x(k) = d(k) + b(k)*xn
      end do
      return
      end

