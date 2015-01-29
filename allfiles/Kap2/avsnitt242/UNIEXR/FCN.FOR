      subroutine fcn(neq,x,y,dy)
      implicit double precision (a-h,o-z)
      dimension y(3),dy(3)
      common /param/ ivelg

c     the blasius equation for two uniform mixing layers
c     coordinate is eta in the upper layer and ksi
c     in the lower layer.
c     ivelg = 1 : upper layer, i.e : eta >= 0
c     ivelg = 2 : lower layer, i.e : ksi >= 0

      dy(1) = y(2)
      dy(2) = y(3)
      dy(3) = -y(1)*y(3) 
      if (ivelg .eq. 2 ) dy(3) = - dy(3)

      return
      end
