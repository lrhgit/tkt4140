% Program ADVECEX1
% Eksplisitt skjema for adveksjonsligningen
% Genererer en av figurene i fig. 6.3 i manus
% Leser inn Courant-tallet C og antall tidskritt nitr
% Bruker dx = 0.01
% Eksempel på brukbare verdipar :
% C = 1 , nitr = 10 
% C = 0.5 , nitr = 20
% C = 0.25 , nitr = 40 
n = 101; % Antall punkt
dx = 1/(n-1);
u = zeros(n,1); u0 = u; x = u;
for k = 1:n
    x(k) = (k - 1)*dx; % 0.0 <= x <= 1.0 
end
% Startverdier for u.
% u = 1.0 for x <= 0.5 ,ellers u = 0
kmax = round(n/2);
for k = 1:kmax
    u(k) = 1.0;  
end
nitr = input('Antall iterasjoner = ?');
C = input('Courant-tall = ?');
xfront = nitr*C*dx + 0.5; % Front of exact solution
for k = 1:n
    if x(k)< xfront   
      u0(k) = 1.0;
  end   
end
for itr = 1: nitr
    a = u(1); b = u(2);
    for j = 2: n-1
        u(j) = b - C*(u(j+1) - a)*0.5;
        a = b; b = u(j+1); 
    end
end
FS = 'FontSize';
clf
plot(x,u,'k',x,u0,'k:');
xlabel('x',FS,14);
ylabel('u(x,t)',FS,14);
st = sprintf('Antall tidskritt = %3.0f',nitr);
title(st,FS,14);
%subplot(312),plot(x,u,'k',x,u0,'k:');
%subplot(313),plot(x,u,'k',x,u0,'k:');
%ylim([-2.0 2.0]);