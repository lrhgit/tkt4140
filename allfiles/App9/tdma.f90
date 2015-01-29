!************************ T D M A *******************************
!                                                               *
!     Solution of a linear system of algebraic equations with   *
!     a tridiagonal matrix of coefficients.(No pivoting)        *
!     Equation no. i :                                          *
!         a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i),         *
!                                            i = 1,2,...n       *
!                       === USE ===                             *
!                                                               *
!                  call tdma(n,a,b,c,d,x)                       *
!                            or                                 *
!                  call tdma(n,a,b,c,d,d)                       *
!                                                               *
!      In the last case, vector d contains the solution.        *
!                                                               *
!                     === INPUT ===                             *
!                                                               *
!     n ....... integer     . Number of equations               *
!     a(1:n) .. real vector . Lower diagonal.Element A(1)       *
!                                            is not used.       *
!     b(1:n) .. real vector . Main diagonal                     *
!     c(1:n) .. real vector . Upper diagonal.Element C(N)       *
!                                            is not used.       *
!     d(1:n) .. real vector . Right hand side of the system.    *
!                                                               *
!                     === OUTPUT ===                            *
!                                                               *
!     x(1:n) .. real vector . The solution vector               *
!                                                               *
!********************** fortran 90 ******************************

      subroutine tdma(n,a,b,c,d,x)
	  implicit none
      integer, intent(in) :: n
      real, intent(in) :: a(n), c(n)
      real, intent(inout), dimension(n) :: b, d
	  real, intent(out) :: x(n)
	  !  --- Local variables ---
	  integer :: i
	  real :: q
      !  --- Elimination ---
      do i = 2,n
         q = a(i)/b(i - 1)
         b(i) = b(i) - c(i - 1)*q
         d(i) = d(i) - d(i - 1)*q
      end do
      ! --- Backsubstitution ---
      q = d(n)/b(n)
      x(n) = q
      do i = n - 1,1,-1
         q = (d(i) - c(i)*q)/b(i)
         x(i) = q
      end do
      return
      end
