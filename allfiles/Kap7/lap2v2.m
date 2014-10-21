% program lap2v2
% Solves the temperature-problem in figure 7.6, section 7.2.2
% Direct solver.
% The temperatures are stored in the vector T.
% This version is an "automatic" version of lap2v1
%
clear
h = 0.25;
nl = 1/h; % No. of unknowns per line
neq = nl^2 ; % Total no. of equations
d0 = zeros(neq,1); % diagonal
d = ones(neq,1); % diagonal
b = d0; % right hand side
% --- Modify b
for k = 0 : nl -1
    b(neq - k) =  - 1;
end
% --- Generate A-matrix ---
for k = -nl:nl
    diagvec(k + nl + 1) = k;
end
A = spdiags([d d0 d0 d -4*d d d0 d0 d],diagvec, neq,neq);
% === Update A ===
for k = 1: nl % superdiagonal no. nl
    A(k,nl + k) = 2;
end
k = 1 ;  % superdiagonal no. 1
while k < neq
    A(k,k + 1) = 2.0 ;
    k = k + nl ;
end
% --- super- and subdiagonal no. 1
k = nl ;
while k < neq
    A(k,k + 1) = 0.0 ;
    A(k + 1,k) = 0.0 ;
    k = k + nl ;
end
% --- solve system ---
T = A\b;