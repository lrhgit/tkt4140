      subroutine fvalue(n,x,fvec,iflag)
c
c     x(1) = r, x(2) = s
c     fvec(1) = psi, fvec(2) = phi
c
      implicit double precision (a-h,o-z)
      double precision ksi,ksimax
      common /share/ h,hmax,eps,nunit,nout,u2u1,etamax,ksimax
      common /param / ivelg
      dimension fe(3),x(n),fk(3),fvec(n)
      external fcn
      data zero,one /0.0d0, 1.0d0/

      neq = 3
      fe(1) = zero
      fe(2) = x(1)
      fe(3) = x(2)
      eta  = zero
      err = eps
      ivelg = 1
      call odex87(neq,fcn,eta,fe,etamax,err,
     &            hmin,hmax,h,nunit,iflag)

      fvec(2) = fe(2) - one
c      write(4,*) ' phi = ',fvec(2)
      fk(1) = zero
      fk(2) = - x(1)
      fk(3) = x(2)
      ksi  = zero
      ivelg = 2

      call odex87(neq,fcn,ksi,fk,ksimax,err,
     &            hmin,hmax,h,nunit,iflag)

      fvec(1) = fk(2) + u2u1
c      write(4,*) ' psi = ',fvec(1)
      return
      end
