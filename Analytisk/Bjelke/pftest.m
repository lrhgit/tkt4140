% program pftest
% Tester funksjonen ellipfk
%
phivec = [5 30 60 85]';
n = length(phivec);
k = 0.5*sqrt(3);
for l = 1:n
    phi = phivec(l)*pi/180;
    [K,F] = ellipfk(phi,k);
    fprintf('F = %20.16e \n',F);
end

