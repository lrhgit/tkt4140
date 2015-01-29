!*************************** B D S O L V *******************************
!                                                                      *
!   Computes the solution x from the linear set of equations           *
!                            a*x = b                                   *
!   where a is a band matrix of coefficients .                         *
!   A triangular factorization with partial pivoting is applied to a.  *
!   Use banfac and bansol instead if you want to solve the system      *
!   with the same a but with different right hand sides b.             *
!                                                                      *
!   bdsolv is a combined fortran version of the algol subroutines      *
!   bandet1 and bansol given in the reference below.                   *
!                                                                      *
!   =======================  Use  ==================================   *
!                                                                      *
!              call bdsolv(n,m1,m2,a,b,x,ifail)                        *
!                                                                      *
!       In order to save space, you may store x in b :                 *
!                                                                      *
!              call bdsolv(n,m1,m2,a,b,b,ifail)                        *
!                                                                      *
!   ======================= Input ==================================   *
!                                                                      *
!   n  - integer    - : Order of a (No. of equations)                  *
!   m1 - integer    - : Number of subdiagonals                         *
!   m2 - integer    - : Number of superdiagonals                       *
!   a - real array  - : Matrix of coefficients dimensioned             *
!                       (1:n,1:m1+m2+1) ; i.e: n*(m1+m2+1) elements.   *
!                       The original contents of a is destroyed        *
!                       in the elimination.                            *
!                       For details on storage, see below.             *
!                                                                      *
!   b - real array - : Vector dimensioned (1:n). Right hand side of    *
!                      the system of equations. Destroyed during the   *
!                      forward substitution.                           *
!                                                                      *
!   ======================= Output  ================================   *
!                                                                      *
!   x - real array - : Vector dimensioned (1:n). The solution vector.  *
!                      If you use the second version of the call       *
!                      given above, the solution is output in the      *
!                      vector b.                                       *
!   ifail -integer - : ifail =  0 -> successful factorization          *
!                      ifail = -1 -> Singular matrix A (Zero pivot)    *
!                                                                      *
!   ===================== Details on storage =======================   *
!                                                                      *
!   The main diagonal is stored in a(i,m1+1),i = 1,2..n.               *
!   The subdiagonal most distant from the maindiagonal is stored in    *  
!   a(i,1), i = 1,2,..n, with a(i,1) = 0 for i = 1,..m1.               *       
!   The next subdiagonal is stored in a(i,2), i = 1,2,..n, with        *
!   a(i,2) = 0 for i = 1,..m1 - 1; etc.                                *
!   The first superdiagonal is stored in a(i,m1+2), i = 1,2,..n,       *
!   with a(i,m1+2) = 0 for i = n-1,n and the last superdiagonal        *              
!   is stored in a(i,m1+m2+1),i = 1,..n, with a(i,m1+m2+1) = 0         *
!   for i = n - m2,..n.                                                *
!                                                                      *
!   ====================  Reference  ==============================    *
!                                                                      *
!   Wilkinson & Reinsch ; "Handbook for Automatic Computation,         *
!                          Vol. II : Linear Algebra",                  *
!                          Springer-Verlag 1971                        *
!                                                                      *
!**************************** fortran 90 *******************************

    subroutine bdsolv(n,m1,m2,a,b,x,ifail)
	implicit none
    integer, intent(in) :: n,m1,m2
	integer, intent(out)  :: ifail
    real, intent(inout) ::  a(n,1:m1+m2+1),b(n)
	real, intent(out) :: x(n)
	!  --- Local variables ---
	integer :: m3,l,i,j,k,i2 
	real :: q, q1
    real, parameter :: zero = 0.0e0

    ifail = 0
    !  --- Reordering ---
    m3 = m1 + 1
    l = m1
    do i = 1,m1
      do j = 1 - i + m3,m2 + m3
         a(i,j - l) = a(i,j)
      end do
      l = l - 1
      do j = m2 - l,m2
         a(i,j+m3) = zero
      end do
    end do
    !  --- Elimination ---
    l = m1
    do k = 1,n
       q = a(k,1)
       i = k
       if(l < n) l = l + 1
       do j = k + 1,l
          q1 = a(j,1)
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
          q = b(k)
          b(k) = b(i)
          b(i) = q
          do j =1,m2 + m3
             q = a(k,j)
             a(k,j) = a(i,j)
             a(i,j) = q
          end do
       endif
       do i = k + 1,l
          q = a(i,1)/a(k,1)
          b(i) = b(i) - q*b(k)
          do j = 2,m2 + m3
             a(i,j-1) = a(i,j) - q*a(k,j)
          end do
		  a(i,m2+m3) = zero
       end do
    end do
    !  --- Backsubstitution ---
    l = - m1
    do i = n,1,-1
       q = b(i)
       x(i) = q
       i2 = i + m1
       do k = 1 - m1,l
          q = q - a(i,k+m3)*x(k + i2)
       end do
       x(i) = q/a(i,1)
       if(l < m2) l = l + 1
    end do
    return
    end
