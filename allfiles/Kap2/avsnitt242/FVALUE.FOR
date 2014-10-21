      subroutine fvalue(phi,psi,s,r)
      implicit double precision (a-h,o-z)
      double precision ksi,ksimax
      common /share/ h,hmax,eps,nunit,nout,u2u1,etamax,ksimax
      common /param / ivelg
      dimension fe(3),fk(3)
      external fcn
      data zero,one /0.0d0, 1.0d0/

      neq = 3
      fe(1) = zero
      fe(2) = r
      fe(3) = s
      eta  = zero
      err = eps
      ivelg = 1
      call odex87(neq,fcn,eta,fe,etamax,err,
     &            hmin,hmax,h,nunit,iflag)

      phi = fe(2) - one

      fk(1) = zero
      fk(2) = - r
      fk(3) = s
      ksi  = zero
      ivelg = 2

      call odex87(neq,fcn,ksi,fk,ksimax,err,
     &            hmin,hmax,h,nunit,iflag)

      psi = fk(2) + u2u1

      return
      end
