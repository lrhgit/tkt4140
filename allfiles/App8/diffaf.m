% ******************** Program diffaf *********************
% This program computes finite difference expressions
% based on Taylor expansions about point no. j which is
% denoted point 0 in the program.
% The programs reads and writes to the screen.
% Examples of use are given in the compendium.
% Reference : S. Nakamura : " Applied Numerical Methods
%                             with Software ",
%                             Prentice-Hall 1991
%*********************************************************

fprintf(' === Difference Approximation Finder ===\n\n');
while true   
   km = input(' Number of points ?');
   if km >= 2 , break, end    
   disp ('   km < 2 . Repeat input!');        
end
% --- Initialize ---
km1 = km +1; km2 = km + 2;
a = zeros(km2,km1);
el = zeros(km,1); x = el; cf = el;
c = zeros(km2,1); d = c;
eps1 = 1.0e-8; eps2 = 1.0e-12;

for k = 1: km
   fprintf('   Index of point no. %2.0f ?',k);
   el(k) = input('');
end
disp(' ');
while true
   kdr = input('  Order of derivative ?');
   if km >= (kdr + 1), break, end   
   disp('   Order of derivative too high');           
end
for k = 1:km2
   for l = 1:km 
      if k == 1
         a(k,l) = 1;         
      else      
         a(k,l) = el(l)^(k - 1);             
      end
   end         
end
z = factorial(kdr);
for k = 1 : km   
   if k == kdr + 1, d(k) = z ; end
end
disp(' ');
x = a(1:km,1:km)\d(1:km);
d(1:km) = x;
c = a(1:km2,1:km)*x;
f = 1000;
for k = 1: km
   u = abs(d(k));
   if u > eps1
      if u < f, f = u ; end    
   end
end
cf = d(1:km)/f;
fprintf('     Difference Scheme \n\n');
s = '+[ %10.5f /( %8.5f H**%1.0f )] Y(%6.3fH)\n';
for k = 1:km
   finv = 1/f;   
   fprintf(s,cf(k),finv,kdr,el(k));      
end
disp(' ');
fprintf('     Error Term \n\n');
for k = 1: km2
   if(abs(c(k)) < eps2), c(k) = 0; end 
end
dd = factorial(km);
s1 = '    ( %15.3f/%10.3f)H**%1.0f  Y^(%1.0f)\n';
for k = km1: km2
   cm = -c(k);
   kmi = k -1;
   nh = kmi - kdr;
   if (k == km1) & (cm ~= 0)
      fprintf(s1,cm,dd,nh,kmi);
   end
   if  k == km2
      fprintf(s1,cm,dd,nh,kmi);
   end
   dd = dd*k;
end

   