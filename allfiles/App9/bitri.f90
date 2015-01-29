!************************* B I T R I  **************************
!                                                               *
!     Solution of a linear system of algebraic equations with   *
!     a bi-tridiagonal matrix of coefficients (no pivoting)     *
!     using the Douglas algorithm. This version does not        *
!     destroy the contents of the given vectors.                *
!                                                               *
!     Equation no. k (two equations):                           *
!                                                               *
!         a1(k)*u(k-1) + a2(k)*v(k-1) + b1(k)*u(k)              *
!       + b2(k)*v(k)   + c1(k)*u(k+1) + c2(k)*v(k+1) = d1(k)    *
!                                                               *
!         a3(k)*u(k-1) + a4(k)*v(k-1) + b3(k)*u(k)              *
!       + b4(k)*v(k)   + c3(k)*u(k+1) + c4(k)*v(k+1) = d2(k)    *
!                                                               *
!                                            k = 1,2,...n       *
!                                                               *
!     This system is equivalent to a block tridiagonal system   *
!     with 2x2 blocks                                           *
!                                                               *
!     Input:                                                    *
!                                                               *
!     n ......... integer  . Number of equations                *
!                                                               *
!     a1 - a4 ... real  vectors dimensioned 1:n                 *
!                 Elements am(1),m = 1,..,4 are not used.       *
!     b1 - b4 ... real vectors dimensioned 1:n                  *
!                                                               *
!     c1 - c4 ... real vectors dimensioned 1:n                  *
!                 Elements cm(n),m = 1,..,4 are not used.       *
!     d1 - d2 ... real vectors dimensioned 1:n                  *
!                 Right hand sides of the system.               *
!                                                               *
!     Output:                                                   *
!                                                               *
!     u and v ... real vectors dimensioned 1:n                  *
!                 The solution of the system                    *
!                                                               *
!************************ for90 *********************************

      subroutine bitri(n,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4, &
                     d1,d2,u,v)
      implicit none
      integer, intent(in) :: n
      real, intent(in), dimension(n) :: a1,a2,a3,a4,b1,b2,b3,b4, &
                                        c1,c2,c3,c4,d1,d2
      real, intent(out) :: u(n), v(n)
	  !   --- Local variables ---
	  integer :: k
      real, dimension(n) :: ba1,ba2,ba3,ba4,ga1,ga2,ga3, &
                            ga4,g1,g2,da1,da2
      real :: my
      real, parameter :: one = 1.0e0

      ba1(1) = b1(1)
      ba2(1) = b2(1)
      ba3(1) = b3(1)
      ba4(1) = b4(1)
      da1(1) = d1(1)
      da2(1) = d2(1)
      my = one/(ba1(1)*ba4(1) - ba2(1)*ba3(1))
      ga1(1) = (ba4(1)*c1(1)  - ba2(1)*c3(1))*my
      ga2(1) = (ba4(1)*c2(1)  - ba2(1)*c4(1))*my
      ga3(1) = (ba1(1)*c3(1)  - ba3(1)*c1(1))*my
      ga4(1) = (ba1(1)*c4(1)  - ba3(1)*c2(1))*my
      g1(1)  = (ba4(1)*da1(1) - ba2(1)*da2(1))*my
      g2(1)  = (ba1(1)*da2(1) - ba3(1)*da1(1))*my
      !  --- Elimination ---
      do k = 2,n
         ba1(k) = b1(k) - a1(k)*ga1(k-1) - a2(k)*ga3(k-1)
         ba2(k) = b2(k) - a1(k)*ga2(k-1) - a2(k)*ga4(k-1)
         ba3(k) = b3(k) - a3(k)*ga1(k-1) - a4(k)*ga3(k-1)
         ba4(k) = b4(k) - a3(k)*ga2(k-1) - a4(k)*ga4(k-1)
         da1(k) = d1(k) - a1(k)*g1(k-1)  - a2(k)*g2(k-1)
         da2(k) = d2(k) - a3(k)*g1(k-1)  - a4(k)*g2(k-1)
         my = one/(ba1(k)*ba4(k) - ba2(k)*ba3(k))
         ga1(k) = (ba4(k)*c1(k) - ba2(k)*c3(k))*my
         ga2(k) = (ba4(k)*c2(k) - ba2(k)*c4(k))*my
         ga3(k) = (ba1(k)*c3(k) - ba3(k)*c1(k))*my
         ga4(k) = (ba1(k)*c4(k) - ba3(k)*c2(k))*my
         g1(k)  =  (ba4(k)*da1(k) - ba2(k)*da2(k))*my
         g2(k)  =  (ba1(k)*da2(k) - ba3(k)*da1(k))*my
      end do
      ! --- Backsubstitution ---
      u(n) = g1(n)
      v(n) = g2(n)
      do k = n-1,1,-1
         u(k) = g1(k) - ga1(k)*u(k+1) - ga2(k)*v(k+1)
         v(k) = g2(k) - ga3(k)*u(k+1) - ga4(k)*v(k+1)
      end do
      return
      end


