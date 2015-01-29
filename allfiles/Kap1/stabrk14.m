% Program stabrk14      
% ~~~~~~~~~~~~~~~~~~
% This program plots the boundary of the region
% of the complex plane governing the maximum 
% step size which may be used for stability of 
% a Runge-Kutta integrator of  order 1 to 4.
%
% npts  - a value determining the number of 
%         points computed on the stability 
%         boundary of an explicit Runge-Kutta 
%         integrator.
% xrang - controls the square window within 
%         which the diagram is drawn. 
%         [ -3, 3, -3, 3] is appropriate for 
%         the fourth order integrator.
%
% User m functions required: none

hold off; close;

npts = 150;
tegnsize = 4;
nordrmax = 4;
tegnsize = 5;
hold on;
for nordr = 1:nordrmax
    r=zeros(npts,nordr); v=1./gamma(nordr+1:-1:2);
    d=2*pi/(npts-1); i=sqrt(-1);
    % Generate polynomial roots to define the 
    % stability boundary
    for j=1:npts
         % polynomial coefficients
        v(nordr+1)=1-exp(i*(j-1)*d); 
        % complex roots
        t=roots(v); r(j,:)=t(:).';
    end
    % Plot the boundary
    rel = real(r(:)); img=imag(r(:)); 
    w = 1.1*max(abs([rel;img]));
    zoom on; plot(rel,img,'.','markersize',tegnsize); 
end     
hold off;
axis([-3,3,-3,3]); axis('square');
xlabel('real-del av \lambdah','FontSize',14);
ylabel('imaginær-del av \lambdah','FontSize',14);
title('Stabilitetsområde for eksplisitte RK-skjema','FontSize',14)
grid on; 
shg
