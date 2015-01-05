function [u,v] = douglas(a2,a3,b1,b2,b3,b4,c2,c3,d1,d2)
%
% Spesialversjon av bitri for bruk i programmet bladoug.
% Da vektorene a1, a4, c1 og c4 = 0 i bladoug, har vi fjernet
% disse vektorene i kallet samt internt.
%
n = length(b1);
u = zeros(size(b1)); v = u;
ga1 = u; ga2 = u; ga3 = u; ga4 = u; da1 = u; da2 = u;
ba1 = u; ba2 = u; ba3 = u; ba4 = u; g1 = u; g2 = u;

ba1(1) = b1(1); ba2(1) = b2(1); ba3(1) = b3(1); ba4(1) = b4(1);
da1(1) = d1(1); da2(1) = d2(1);

my = 1.0/(ba1(1)*ba4(1) - ba2(1)*ba3(1));
ga1(1) = - ba2(1)*c3(1)*my;
ga2(1) = ba4(1)*c2(1)*my;
ga3(1) = ba1(1)*c3(1)*my;
ga4(1) = - ba3(1)*c2(1)*my;
g1(1)  = (ba4(1)*da1(1) - ba2(1)*da2(1))*my;
g2(1)  = (ba1(1)*da2(1) - ba3(1)*da1(1))*my;

%     ===== Elimination ====

for k = 2:n
   ba1(k) = b1(k) - a2(k)*ga3(k-1);
   ba2(k) = b2(k) - a2(k)*ga4(k-1);
   ba3(k) = b3(k) - a3(k)*ga1(k-1);   
   ba4(k) = b4(k) - a3(k)*ga2(k-1);
   da1(k) = d1(k) - a2(k)*g2(k-1);
   da2(k) = d2(k) - a3(k)*g1(k-1);
   my = 1.0/(ba1(k)*ba4(k) - ba2(k)*ba3(k));
   ga1(k) = -ba2(k)*c3(k)*my;
   ga2(k) =  ba4(k)*c2(k)*my;
   ga3(k) =  ba1(k)*c3(k)*my;
   ga4(k) = -ba3(k)*c2(k)*my;
   g1(k) =  (ba4(k)*da1(k) - ba2(k)*da2(k))*my;
   g2(k) =  (ba1(k)*da2(k) - ba3(k)*da1(k))*my;
end

%    ==== Backsubstitution ====

u(n) = g1(n); v(n) = g2(n);
for k = n-1:-1:1
   u(k) = g1(k) - ga1(k)*u(k+1) - ga2(k)*v(k+1);
   v(k) = g2(k) - ga3(k)*u(k+1) - ga4(k)*v(k+1);
end 

      