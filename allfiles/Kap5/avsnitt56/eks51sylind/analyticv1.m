% program analyticv1
% Beregner den analytiske løsningen av diff. ligningen
%   dw/dt = w''(r) + w'(r)/r , w = w(r,t), 0 < r < 1
dr = 0.1; % Skrittlengde
t = 0; % Tid
rvec = (dr : dr: 1-dr)';
N = length(rvec);
w = zeros(N,1);
for k = 1:N
    r = rvec(k);
    % Beregner 1. ledd forseg selv
    n = 1;
    lam1 = j0zero(n);
    arg1 = r*lam1;
    term1 = exp(-t*lam1^2)*besselj(0,arg1)/(besselj(1,lam1)*lam1^3);  
    sum1 = term1; sumtol = 1.0e-5; test = 1;
    while test > sumtol
        n = n + 1;
        lamn = j0zero(n);
        arg = r*lamn;
        term = exp(-t*lam1^2)*besselj(0,arg)/(besselj(1,lamn)*lamn^3);       
        sum1 = sum1 + term;
        test = abs(term/term1);
    end
    w(k) = 8*sum1;
end
rvec = [0;rvec;1];
w = [1;w;0];
[rvec w]