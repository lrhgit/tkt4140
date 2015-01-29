function x = tridiag(a,b,c,d)
% Solution of a tri-diagonal matrix 
% using the Thomas-algorithm.
% Column-wise numeration of the coefficients:
% a(i-1)*x(i-1) + b(i)*x(i) + c(i+1)*x(i+1) = d(i)
% No. of equations given by the length
% of the main diagonal b
n = length(b);
x = zeros(size(b));
% === Elimination ===
for k = 2:n
   q = a(k-1)/b(k-1);
   b(k) = b(k) - c(k)*q;
   d(k) = d(k) - d(k-1)*q;
end
% === Backsubstitution ===
x(n) = d(n)/b(n);
for k = n-1 :-1 :1
    x(k) = (d(k) - c(k+1)*x(k+1))/b(k);
end
