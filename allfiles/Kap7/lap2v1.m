% program lap2v1
% Solves case in figure 7.6, section 7.2.2
% Direct solver
clear
n = 16;
h = 0.25; h2 = h*h;
d0 = zeros(n,1); % diagonal
d = ones(n,1); % diagonal
b = d0; % right hand side
% --- Modify b
for k = 13:16
    b(k) =  - 1;
end
% --- generate A-matrix ---
A = spdiags([d d0 d0 d -4*d d d0 d0 d],[-4 -3 -2  -1 0 1 2 3 4], n,n);
% === Update A ===
A(1,2) = 2; A(1,5) = 2; A(2,6) = 2; A(3,7) = 2;
A(4,5) = 0; A(4,8) = 2; A(5,4) = 0; A(5,6) = 2; 
A(8,9) = 0; A(9,8) = 0; A(9,10) = 2;
A(12,13) = 0; A(13,12) = 0; A(13,14) = 2; 

% --- solve system ---
T = A\b;