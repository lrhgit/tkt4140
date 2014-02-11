% program lap2v3
% Solves the temperature-problem in figure 7.6, section 7.2.2
% Direct solver.
% The temperatures are stored in the vector T.
% In this version we can select 4 different nets
% by specifying the parameter net.
%
% neq = total number of unknowns = no. of equations
% net = 1 -> h = 0.25             -> neq = 16
% net = 2 -> h = 0.25/2 = 0.125   -> neq = 64
% net = 3 -> h = 0.25/4 = 0.0625  -> neq = 256
% net = 4 -> h = 0.25/8 = 0.03125 -> neq = 1024
%
clear
h = 0.25;
net = 3; % Select net
n = net - 1;
hn = h/(2^n);
nl = 1/hn; % No. of unknowns per line
neq = nl^2 ; % Total no. of equations
fprintf('Gridsize h =  %12.5e \n',hn);
fprintf('No. of equations =  %10.0f \n',neq);
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
switch net
    case 1
        A = spdiags([d d0 d0 d -4*d d d0 d0 d],diagvec, neq,neq);
    case 2
        A = spdiags([d d0 d0 d0 d0 d0 d0 d -4*d d d0 d0 d0 d0 d0 d0 d],diagvec, neq,neq);
    case 3
        A = spdiags([d d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d -4*d ... 
        d d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d],diagvec, neq,neq);
    case 4
        A = spdiags([d d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 ...
                     d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d -4*d ... 
                     d d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 ...
                     d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d],diagvec, neq,neq);
    otherwise
        error('Parameter net out of range')
end
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