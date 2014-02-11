function zj0  = j0zero3(s)
% Function j0zero computes root number s, s = 1,2,.., 
% of the Besselfunction J0 where zj0 is the root,
% using the Newton-Raphsons method.
% Relative error ca. 1.0e-15

persistent root;
if isempty(root)
    % --- Table for the first five zeros
    root(1) = 2.404825557695773;
    root(2) = 5.520078110286311;
    root(3) = 8.653727912911012;
    root(4) = 11.79153443901428;
    root(5) = 14.93091770848779;
end
if (s <= 5) & (s >= 1)
    zj0 = root(s);
    return
end
% --- Compute an initial value of the root
b0 = (s - 0.25)*pi;
b08 = 0.125/b0; b082 = b08^2; 
z0 = b0 + b08*(1 - b082*(124/3 - b082*120928/15)); 
if (s > 30)
    zj0 = z0;
    return
end
% ---  5 < s < 30 ---
    dz0 = besselj(0,z0)/besselj(1,z0);
    z0 = z0 + dz0;
    zj0 = z0;
