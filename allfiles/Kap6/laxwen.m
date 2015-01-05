% Program LAXWEN
% Lax-Wendroff skjema for adveksjonsligningen
% for gitt verdi av C.
% Bruker dx = 0.01
% Eksempel på verdier :
% C = 1 , nitr = 10 
% C = 0.5 , nitr = 20
% C = 0.25 , nitr = 40 
clear; close;
n = 101; % Antall punkt
dx = 1/(n-1);
u = zeros(n,1); x = u;
u0 = u;
for k = 1:n
    x(k) = (k - 1)*dx; % 0.0 <= x <= 1.0 
end
xfront =  0.6; % Front of exact solution
for k = 1:n
    if x(k)< xfront   
      u0(k) = 1.0;      
   end
end
kmax = round(n/2);
ylim([-0.5 1.5]); 
C = 1;   
u = zeros(n,1);
for k = 1:kmax       
    u(k) = 1.0;        
end
nitr = round(50/C);
k1 = 1 + C; k2 = 1 - C;
for itr = 1: nitr- 1   
    a = u(1); b = u(2);        
    for j = 2: n-1           
        u(j) = k1*k2*b + C*(-k2*u(j+1) + k1*a)*0.5;                
        a = b; b = u(j+1);                    
    end    
end
FS = 'FontSize';
LW = 'LineWidth';
plot(x,u,'k',x,u0,'k:',LW,1);
ylim([-0.5 1.5]); 
%grid
% plot(x,u,'k',x,u0,'k:');
% ylim([-0.25 1.25]); 
%plot(x,u,'k',x1,y1,'k:',x2,y2,'k:',x3,y3,'k:',LW,1.5);
ylabel('u',FS,14,'Rotation',0);   
xlabel('x',FS,14);        
st = sprintf('C = %4.2f',C);
title(st,FS,12);   

