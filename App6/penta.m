function x = penta(e,a,b,c,f,d)
%                                                              
%     Solution of a linear system of algebraic equations with   
%     a pentadiagonal matrix of coefficients.(No pivoting)      
%     Equation no. i :                    
%     e(i)*x(i-2)+ a(i)*x(i-1) + b(i)*x(i) +                    
%     c(i)*x(i+1) + f(i)*x(i+2) = d(i), i = 1,2,...n            
%            where n is the number of equations                       
%
%                     === Input === 
%     Vectors may be in either row- or column-form
%
%     e(1:n) .. Lowest diagonal. e(1 and e(2) are not used          
%     a(1:n) .. Lower diagonal. a(1)is not used.     
%     b(1:n) .. Main diagonal                      
%     c(1:n) .. Upper diagonal. c(n)is not used.            
%     f(1:n) .. Uppermost diagonal. f(n-1)and f(n) are not used.       
%     d(1:n) .. Right hand side of the system. 
%
%                     === Output ===  
%     x(1:n) .. The solution vector                
%                                                                     
n = length(b);
x = zeros(size(b));
% === Elimination
for k = 2 : n-1 
   q = a(k)/b(k-1);
   b(k) = b(k) - q*c(k-1);
   c(k) = c(k) - q*f(k-1);
   d(k) = d(k) - q*d(k-1);
   q = e(k+1)/b(k-1);
   a(k+1) = a(k+1) - q*c(k-1);
   b(k+1) = b(k+1) - q*f(k-1);
   d(k+1) = d(k+1) - q*d(k-1);
end
% === Backsubstitution
q = a(n)/b(n-1);
b(n) = b(n) - q*c(n-1);
x(n) = (d(n) - q*d(n-1))/b(n);
x(n-1) = (d(n-1) - c(n-1)*x(n))/b(n-1);
for k = n-2 : -1 :1   
   x(k) = (d(k) - f(k)*x(k+2) - c(k)*x(k+1))/b(k); 
end
