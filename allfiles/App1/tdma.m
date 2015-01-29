function x = tdma(a,b,c,d)
% Solution of a linear system of algebraic 
% equations with a tri-diagonal matrix of coefficients
% using the Thomas-algorithm (No pivoting).
% Equation no. i :
%  a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i), i = 1,..,n
%  where n is the number of equations.
%
%         === Input ===
% Vectors may be in either row - or column-form
%
%  a(1:n) .. Lower diagonal. a(1) is not used
%  b(1:n) .. Main diagonal. 
%  c(1:n) .. Upper diagonal. c(n) is not used
%  d(1:n) .. Right hand side of the system.
%
%         === Output ===
%  x(1:n) .. The solution vector

n = length(b);
x = zeros(size(b));
% === Elimination ===
for k = 2:n
   q = a(k)/b(k-1);
   b(k) = b(k) - c(k-1)*q;
   d(k) = d(k) - d(k-1)*q;
end
% === Backsubstitution ===
q = d(n)/b(n);
x(n) = q;
for k = n-1 :-1 :1
   q = (d(k) - c(k)*q)/b(k);
   x(k) = q;
end
