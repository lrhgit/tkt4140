      function dy = fcn(x,y)
%       dimension y(3),dy(3)
%       common /param/ ivelg

%     The Blasius equation for two uniform mixing layers
%     The coordinate is eta in the upper layer and ksi
%     in the lower layer.
%     ivelg = 1 : upper layer, i.e : eta >= 0
%     ivelg = 2 : lower layer, i.e : ksi >= 0

      dy(1) = y(2);
      dy(2) = y(3);
      dy(3) = -y(1)*y(3) 
      if (ivelg .eq. 2 )
          dy(3) = - dy(3);
      end