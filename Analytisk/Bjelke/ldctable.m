% program ldctable
% Large deflection of a cantilever beam.
% For given loadfactors alpha the program computes
% the end-deflection delta and the end-slope theta0 
% together with the horisontal length lh.
% The program calls the function kalpha.
% theta0 in degrees.
clear
alphavec = (0:0.5:5)';
n = length(alphavec);
fprintf('alpha   theta0   delta      lh \n')
for l = 1:n
    alpha = alphavec(l);
    k = kalpha(alpha);
    t0 = asin(2*k^2 - 1);
    theta0 = t0*180/pi;
    u1 = asin(1/(sqrt(2)*k));
    [E,Eu1] = ellipek(u1,k);
    if alpha > 1.0e-8    
        delta = 1 - 2*(E - Eu1)/alpha;   
        lh = sqrt(2*sin(t0))/alpha;
    else
        delta = 0;
        lh = 0;
        theta0 = 0;
    end
    fprintf(' %4.1f %8.3f %8.5f %8.5f \n',alpha,theta0,delta,lh);
end