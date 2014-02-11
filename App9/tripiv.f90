!************************** T R I P I V **************************
!                                                                *
!     Solution of a linear system of algebraic equations with    *
!     a tridiagonal matrix of coefficients using pivoting.       *
!     The input-vectors are destroyed in the elimination.        *
!                                                                *
!     Equation no. i :                                           *
!                                                                *
!         a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i),          *
!                                            i = 1,2,...n        *
!                                                                *
!  ========================= USE =============================   *
!                                                                *
!             call tripiv(n,a,b,c,d,x,ifail)                     *
!                                                                *
!  ======================== INPUT ============================   *
!                                                                *
!  n ....... integer     . Number of equations                   *
!  a(1:n) .. real vector . Lower diagonal. Element a(1)          *
!                                          is not used           *
!  b(1:n) .. real vector . Main diagonal                         *
!  c(1:n) .. real vector . Upper diagonal. Element c(n)          *
!                                          is not used           *
!  d(1:n) .. real vector . Right hand side of the system.        *
!                                                                *
!  ========================= OUTPUT ==========================   *
!                                                                *
!  x(1:n) .. real vector . The solution vector                   *
!  ifail  .. integer     . ifail = 0  -> successful elimination  *
!                          ifail = -1 -> singular matrix         *     
!                                                                *
!**************************** for90 ******************************

      subroutine tripiv(n,a,b,c,d,x,ifail)
	  implicit none
      integer, intent(in) :: n
	  integer, intent(out) :: ifail
      real, intent(inout), dimension(n) :: a,b,c,d
	  real, intent(out) :: x(n)
	  !  --- Local variables ---
      integer :: i,j,k,l
	  real :: q,q1
      real, parameter :: zero = 0.e0

      ifail = 0
      ! --- Reordering ---
      a(1) = b(1)
      b(1) = c(1)
      c(1) = zero
      ! --- Elimination ---
      l = 1
      do k = 1,n
         q = a(k)
         i = k
         if(l < n) l = l + 1
         do j = k + 1,l
            q1 = a(j)
            if(abs(q1) > abs(q)) then
               q = q1
               i = j
            endif
         end do
         if( q == zero) then
            ifail = - 1
            return
         endif
         if(i /= k) then
            q    = d(k)
            d(k) = d(i)
            d(i) = q
            q    = a(k)
            a(k) = a(i)
            a(i) = q
            q    = b(k)
            b(k) = b(i)
            b(i) = q
            q    = c(k)
            c(k) = c(i)
            c(i) = q
         endif
         do i = k + 1,l
            q    = a(i)/a(k)
            d(i) = d(i) - q*d(k)
            a(i) = b(i) - q*b(k)
            b(i) = c(i) - q*c(k)
            c(i) = zero
         end do
      end do
      !  --- Backsubstitution ---
      x(n) = d(n)/a(n)
      x(n-1) = (d(n-1) - b(n-1)*x(n))/a(n-1)
      do i = n-2,1,-1
         q = d(i) - b(i)*x(i+1)
         x(i) = (q - c(i)*x(i+2))/a(i)
      end do
      return
      end
