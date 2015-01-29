!************************* b i t r i c **************************
!                                                               *
!     Solution of a linear system of algebraic equations with   *
!     a bi-tridiagonal matrix of coefficients (no pivoting)     *
!     using the Douglas algorithm. This version uses the space  *
!     of the given vectors as workspace and accordingly         *
!     destroys the original contents of those vectors.          *
!                                                               *
!     Equation no. i (two equations):                           *
!                                                               *
!         a1(i)*u(i-1) + a2(i)*v(i-1) + b1(i)*u(i)              *
!       + b2(i)*v(i)   + c1(i)*u(i+1) + c2(i)*v(i+1) = d1(i)    *
!                                                               *
!         a3(i)*u(i-1) + a4(i)*v(i-1) + b3(i)*u(i)              *
!       + b4(i)*v(i)   + c3(i)*u(i+1) + c4(i)*v(i+1) = d2(i)    *
!                                                               *
!                                            i = 1,2,...n       *
!                                                               *
!     This system is equivalent to a block tridiagonal system   *
!     with 2x2 blocks                                           *
!                                                               *
!     Input:                                                    *
!                                                               *
!     n ......... integer  . number of equations                *
!                                                               *
!     a1 - a4 ... real vectors dimensioned 1:n                  *
!                 elements am(1),m = 1,..,4 are not used.       *
!     b1 - b4 ... real vectors dimensioned 1:n                  *
!                                                               *
!     c1 - c4 ... real vectors dimensioned 1:n                  *
!                 elements cm(n),m = 1,..,4 are not used.       *
!     d1 - d2 ... real vectors dimensioned 1:n                  *
!                 right hand sides of the system.               *
!                                                               *
!     Note. the original values of the coefficients             *
!           are destroyed during the process of elimination.    *
!                                                               *
!     Output:                                                   *
!                                                               *
!     u and v ... real vectors dimensioned 1:n                  *
!                 the solution of the system                    *
!                                                               *
!****************************************************************

      subroutine bitric(n,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4, &
                        d1,d2,u,v)
      implicit none  
      integer, intent(in) :: n
      real, intent(inout), dimension(n) :: a1,a2,a3,a4,b1,b2,b3,b4, &
                                           c1,c2,c3,c4,d1,d2
      real, intent(out) :: u(n), v(n)
	  !   --- Local variables ---
	  integer :: i
      real :: my
      real, parameter :: one = 1.0e0

      my = one/(b1(1)*b4(1) - b2(1)*b3(1))
      a1(1) = (b4(1)*c1(1) - b2(1)*c3(1))*my
      a2(1) = (b4(1)*c2(1) - b2(1)*c4(1))*my
      a3(1) = (b1(1)*c3(1) - b3(1)*c1(1))*my
      a4(1) = (b1(1)*c4(1) - b3(1)*c2(1))*my
      u(1) = (b4(1)*d1(1) - b2(1)*d2(1))*my
      v(1) = (b1(1)*d2(1) - b3(1)*d1(1))*my
      ! --- Elimination ---
      do i = 2,n
          b1(i) = b1(i) - a1(i)*a1(i-1) - a2(i)*a3(i-1)
          b2(i) = b2(i) - a1(i)*a2(i-1) - a2(i)*a4(i-1)
          b3(i) = b3(i) - a3(i)*a1(i-1) - a4(i)*a3(i-1)
          b4(i) = b4(i) - a3(i)*a2(i-1) - a4(i)*a4(i-1)

          d1(i) = d1(i) - a1(i)*u(i-1) - a2(i)*v(i-1)
          d2(i) = d2(i) - a3(i)*u(i-1) - a4(i)*v(i-1)

          my = one/(b1(i)*b4(i) - b2(i)*b3(i))
          a1(i) = (b4(i)*c1(i) - b2(i)*c3(i))*my
          a2(i) = (b4(i)*c2(i) - b2(i)*c4(i))*my
          a3(i) = (b1(i)*c3(i) - b3(i)*c1(i))*my
          a4(i) = (b1(i)*c4(i) - b3(i)*c2(i))*my

          u(i) = (b4(i)*d1(i) - b2(i)*d2(i))*my
          v(i) = (b1(i)*d2(i) - b3(i)*d1(i))*my
      end do
      ! --- Backsubstitution ---
      do i = n-1,1,-1
          u(i) = u(i) - a1(i)*u(i+1) - a2(i)*v(i+1)
          v(i) = v(i) - a3(i)*u(i+1) - a4(i)*v(i+1)
      end do
      return
      end
