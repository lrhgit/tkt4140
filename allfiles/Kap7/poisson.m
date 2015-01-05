% program poisson
% Solves example in fig. 7.4 and 7.5, section 7.2.1
clear
n = 15;
h = 0.1; h2 = h*h;
d0 = zeros(n,1); % diagonal
d = ones(n,1); % diagonal
b = -h2*ones(n,1); % right hand side
% --- generate A-matrix ---
A = spdiags([d0 d0 d0 d -4*d d d0 d0 d0],[-4 -3 -2  -1 0 1 2 3 4], n,n);
% === Update A ===
% --- sub-diagonals
A(6,2) = 2; A(7,3) = 1; A(8,4) = 1; A(9,5) = 1;
A(10,7) = 2; A(11,8) = 1; A(12,9) = 1; 
A(13,11) = 2; A(14,12) = 1; 
A(5,4) = 2; A(6,5) = 0; A(9,8) = 2; A(10,9) = 0; A(12,11) = 2;
A(13,12) = 0; A(14,13) = 2; A(15,14) = 4;
% --- super-diagonals ---
A(1,2) = 2; A(5,6) = 0; A(6,7) = 2; A(9,10) = 0;
A(10,11) = 2; A(12,13) = 0; A(13,14) = 2;
A(11,13) = 1; A(12,14) = 1;
A(7,10 ) = 1; A(8,11) = 1; A(9,12) = 1;
A(2,6) = 1; A(3,7) = 1; A(4,8) = 1; A(5,9) = 1;
% --- solve system ---
u = A\b;