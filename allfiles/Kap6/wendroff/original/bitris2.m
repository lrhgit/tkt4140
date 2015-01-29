function [u,v] = bitris2(d1,d2)
global b1 b4 c4;
n = length(d1);
u = zeros(size(d1)); v = u;
ga2 = u; ga4 = u; da1 = u; da2 = u;
ba2 = u; ba4 = u; g1 = u; g2 = u;
ba2(1) = -1; ba4(1) = b4;
da1(1) = d1(1); da2(1) = d2(1);
%
my = 1.0/(b1*b4 + 1);
ga2(1) = -(b4 - c4)*my;
ga4(1) = (b1*c4 + 1)*my;
g1(1)  = (b4*da1(1) + da2(1))*my;
g2(1)  = (b1*da2(1) - da1(1))*my;
%     ===== Elimination ====
for k = 2:n
   ba2(k) = -1 - ga2(k-1);   
   ba4(k) = b4 - ga2(k-1);
   da1(k) = d1(k) - g1(k-1);
   da2(k) = d2(k) - g1(k-1);
   my = 1.0/(b1*b4 - ba2(k));
   ga2(k) = -(ba4(k) +  ba2(k)*c4)*my;
   ga4(k) = (b1*c4 +  1)*my;
   g1(k) =  (ba4(k)*da1(k) - ba2(k)*da2(k))*my;
   g2(k) =  (b1*da2(k) - da1(k))*my;
end
%    ==== Backsubstitution ====
u(n) = g1(n); v(n) = g2(n);
for k = n-1:-1:1
   u(k) = g1(k) - ga2(k)*v(k+1);
   v(k) = g2(k) - ga4(k)*v(k+1);
end 

      