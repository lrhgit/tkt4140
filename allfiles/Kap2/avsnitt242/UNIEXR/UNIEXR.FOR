      program uniexr

c
c     x(1) = r = f'(0) , x(2) = s = f"(0)
c
      parameter (n = 2 ,lw = 19)
      implicit double precision (a-h,o-z)
      double precision ksimax
      dimension x(n), fvec(n), w(lw)
      common /share/ h,hmax,eps,nunit,nout,u2u1,etamax,ksimax
      external fcn, fvalue
      data ninp /3/

c    ===== i n p u t  =====

      iopt = 2
      nprint = 0
      nout = 4
      open(nout,file='uniexr.out')
      open(ninp,file='uniexr.inp')
      write (nout,9990)
      write (nout,9980)
      read (ninp,9810)
      read (ninp,9800) eps,h,hmax,nunit
      if ((nunit .ge. 0) .and. (nunit .ne. nout)) then
         open (nunit,file='solout.out')
      endif
      
      write (nout,9970) eps,h,hmax,nunit
      read (ninp,9780)
      read (ninp,9930) u2u1,etamax,ksimax,tol,x(1),x(2)

      write (nout,9830) u2u1,etamax,ksimax,tol,x(1),x(2)

      write(nout,9960)
c
c     solve nonlinear equations
c
      call dnsqe (fvalue,fcn,iopt,n,x,fvec,tol,nprint,info,w,lw)
      if (info .ne. 1 )write (nout,9770) info
      if (info .ne. -1) then
          write (nout,9760) fvec(1),fvec(2)
          write (nout,9750) x(1) ,x(2)
      endif
   60 stop
 9990 format(1h /
     1       11x,'* * * * * * * * * * * * * * * * * * * * * * * * *'/
     2       11x,'*                                               *'/
     3       11x,'*          uniform mixing of layers             *'/
     4       11x,'*                                               *'/
     5       11x,'*          exact zeros using dnsqe              *'/
     6       11x,'*                                               *'/
     7       11x,'* * * * * * * * * * * * * * * * * * * * * * * * *'/)
 9980 format(/21x,'*** d a t a  i n p u t ***'/)
 9970 format(11x,'max. allowable integration error eps ... ',1pe9.2e2/
     &       11x,'initial grid spacing h ................. ',0pf9.3/
     &       11x,'max. grid spacing hmax ................. ',0pf9.3/
     &       11x,'unit no. for output from solout ........ ',i9/)
 9960 format(/22x,'*** d a t a  o u t p u t ***'/)
 9930 format(5(40x,f15.0/),40x,f15.0)
 9830 format(11x,'velocity ratio u2/u1 ................... ',f9.3/
     &       11x,'max. value of eta etamax ............... ',f9.3/
     &       11x,'max. value of ksi ksimax ................',f9.3/
     &       11x,'max. error in dnsqe tol ................ ',1pe9.2e2/
     &       11x,'guessed value of r = x(1) .............. ',0pf9.3/
     &       11x,'guessed value of s = x(2) .............. ',0pf9.3/)
 9810 format(///////)
 9800 format(3(40x,f15.0/),40x,i15)
 9780 format(/)
 9770 format(15x,' info = ',i4/)
 9760 format(15x,'psi = ',1pe12.4,3x,'phi = ',1pe12.4/)
 9750 format(15x,'r = ',1pe15.7,3x,'s = ',1pe15.7)
      end
