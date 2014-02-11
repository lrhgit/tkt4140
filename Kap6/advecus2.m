% Program ADVECUS2
% Oppstrømskjema skjema for adveksjonsligningen
% Genererer alle figurene i fig. 6.3 i manus
% ved å velge C = 1.0, 0.5 0g 0.25
% Bruker dx = 0.01
% Velger følgende verdipar :
% C = 1 , nitr = 10 
% C = 0.5 , nitr = 20
% C = 0.25 , nitr = 40 
n = 101; % Antall punkt
dx = 1/(n-1);
u = zeros(n,1); x = u;
for k = 1:n
    x(k) = (k - 1)*dx; % 0.0 <= x <= 1.0 
end
x1 = [0.2 0.6]; y1 = [1.0 1.0]; % Brukes for å tegne
x2 = [0.6 0.6]; y2 = [0.0 1.0]; % den analytiske løsningen
x3 = [0.6 0.7]; y3 = [0.0 0.0]; 
xfront =  0.6; % Front of exact solution
% for k = 1:n
%     if x(k)< xfront   
%       u0(k) = 1.0;      
%    end
% end
kmax = round(n/2);
nplot = 0;

for C = [1.0 0.5 0.25]
   nplot = nplot + 1;
   if nplot == 1
       u1 = [0.0 0.6]; v1 = [1.0 1.0]; % Brukes for å tegne
       u2 = [0.6 0.6]; v2 = [0.0 1.0]; % løsningen for C =1
       u3 = [0.6 1.0]; v3 = [0.0 0.0]; % som er den analytiske løsningen
   end      
   % Startverdier for u
   % u = 1.0 for x <= 0.5 ,ellers u = 0
   u = zeros(n,1);
   for k = 1:kmax       
      u(k) = 1.0;        
   end   
   nitr = round(10/C);
   for itr = 1: nitr    
      a = u(1); b = u(2);    
      for j = 2: n-1       
         u(j) = (1 - C)*b + C*a;        
         a = b; b = u(j+1);          
      end      
   end
   FS = 'FontSize'; LW = 'LineWidth';
   subplot(3,1,nplot)
   % plot(x,u,'k',x,u0,'k:');
   if nplot == 1   
       plot(u1,v1,'k',u2,v2,'k',u3,v3,'k',LW,1.5);
   else      
       plot(x,u,'k',x1,y1,'k:',x2,y2,'k:',x3,y3,'k:',LW,1.5);
   end
   ylabel('u',FS,14,'Rotation',0); 
   ylim([-0.5 1.5]); 
    %  xlabel('x',FS,14);     
  
%    st = sprintf('Antall tidskritt = %3.0f',nitr);
%    title(st,FS,14);   
end
