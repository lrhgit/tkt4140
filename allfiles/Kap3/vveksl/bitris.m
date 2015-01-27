function [u,v]= bitris(a1,a3,b1,b2,b3,b4,c2,c4,d1,d2)
      n = length(b1);
      u = zeros(size(b1));
      v = u;

      my = 1/(b1(1)*b4(1) - b2(1)*b3(1));
      a1(1) = (b4(1)*c2(1) - b2(1)*c4(1))*my;
      a3(1) = (b1(1)*c4(1) - b3(1)*c2(1))*my;

      u(1) = (b4(1)*d1(1) - b2(1)*d2(1))*my;
      v(1) = (b1(1)*d2(1) - b3(1)*d1(1))*my;
     % ===== Elimination ====
      for k = 2:n
          b2(k) = b2(k) - a1(k)*a1(k-1); 
          b4(k) = b4(k) - a3(k)*a1(k-1); 

          d1(k) = d1(k) - a1(k)*u(k-1); 
          d2(k) = d2(k) - a3(k)*u(k-1); 

          my = 1/(b1(k)*b4(k) - b2(k)*b3(k));
          a1(k) = (b4(k)*c2(k) - b2(k)*c4(k))*my;
          a3(k) = (b1(k)*c4(k) - b3(k)*c2(k))*my;

          u(k) = (b4(k)*d1(k) - b2(k)*d2(k))*my;
          v(k) = (b1(k)*d2(k) - b3(k)*d1(k))*my;
      end
    % ==== Backsubstitution ====
      for k = n-1:-1:1
          u(k) = u(k) - a1(k)*v(k+1);
          v(k) = v(k) - a3(k)*v(k+1);
      end