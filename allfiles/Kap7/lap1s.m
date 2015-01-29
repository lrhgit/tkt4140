% program lap1s
% Solves example in fig. 7.3, section 7.2.1 and taking
% account of symmetry.
clear
n = 10;
d = ones(n,1); %  diagonal
b = zeros(n,1); % right hand side
%  --- Update b ---
b(5) = -100; b(10) = -100; 
% --- generate A-matrix ---
A = spdiags([2*d d -4*d d d],[-5 -1 0 1 5], n,n);
% --- Update A ---
A(5,6) = 0; A(6,5) = 0; 
% --- solve system ---
T = A\b;