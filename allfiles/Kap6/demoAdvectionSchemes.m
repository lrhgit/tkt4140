% Program demoAdvectionSchemes
% Demo av dei ulike skjema som er demonstrert i kap6 for adveksjonslikninga
% Leif Rune Hellevik

clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

n = 101; % Antall punkt
xmax=1;
dx = xmax/(n-1);
u = zeros(n,1); u0 = u; u2=u; uw=u; ulw=u;
x=linspace(0,xmax,n);

% Startverdiar for u.
% u = 1.0 for x <= 0.5 ,ellers u = 0
kmax = round(n/4);
u(1:kmax) = 1.0;
u2=u;
uw=u;
ulw=u;


% nitr = input('Antall iterasjoner = ?');
% C = input('Courant-tall = ?');
nitr = 160;
C = .35;

umin =-1;
umax = 3;

method='bal';

k1=1+C;
k2=1-C;

for itr = 1: nitr
    xfront = itr*C*dx + kmax*dx; % Front of exact solution
    uptofronti= find(x<=xfront); % indexes of x < xfront
    u0(uptofronti) = 1;          % exact solution

    um  = u(1);  uj  = u(2);
    uwm = uw(1); uwj = uw(2);
    ulwm = ulw(1); ulwj = ulw(2);
  
    
    for j = 2: n-1
        u(j) = uj - C*(u(j+1) -  um)*0.5;
        um = uj; 
        uj = u(j+1);
        
        uw(j) = (1 - C)*uwj + C*uwm;
        uwm = uwj; 
        uwj = uw(j+1);
        
       
        ulw(j) = k1*k2*ulwj + C*(-k2*ulw(j+1) + k1*ulwm)*0.5;                
        ulwm = ulwj; ulwj = ulw(j+1); 
        
        u2(j) = u2(j) - C*(u2(j+1)-u2(j-1))/2;
    end
    
    switch lower(method)
        case{'exact'}
             plot(x,u0);
             legend('Exact');
        
        case{'ftcs'}
            plot(x,u,x,u0,':');
            legend('FTCS','Exact');
        
        case{'upwind'}
            plot(x,uw,x,u0,':');
            legend('Upwind','Exact');

        case{'laxwendroff'}
            plot(x,ulw,x,u0,':');
            legend('LaxWendroff','Exact');
            
        otherwise
            plot(x,uw,x,ulw,x,u0,':');
            legend('Upwind','LaxW','Exact');
    end
    
    set(gca,'ylim',[umin umax]);
    xlabel('x');
    ylabel('u(x,t)');
    st = sprintf('Antall tidskritt = %3.0f',nitr);
    title(st);

    
    pause(0.01);
end

%subplot(312),plot(x,u,'k',x,u0,'k:');
%subplot(313),plot(x,u,'k',x,u0,'k:');
%ylim([-2.0 2.0]);