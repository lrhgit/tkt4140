% Program eksempel3
% Beregning tabell under eksempel 3, del 5.
%
clear
% n = antall ligninger
fprintf('     n        k2       lgk2     epsmax     epsnorm \n\n');
for n = [10 1000 10000 100000 500000]
A = spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],[-1 0 1],n,n);
d = zeros(n,1);
xa = d; 
d(n) = 1;
x = A\d; 
fac = 1/(n+1);
for k = 1:n
    xa(k) = k*fac;% Analytisk løsning
end
f = cos(pi/(n+1));
k2 = (1 + f)/(1 - f); % Tilstandstall
lgk2 = log10(k2);
epsmax = max(abs((xa - x)./xa));
epsnorm = norm(xa - x)/norm(xa);

fprintf('%8.0f %10.2e %7.2f %11.2e %11.2e \n',n,k2,lgk2,epsmax,epsnorm);
end