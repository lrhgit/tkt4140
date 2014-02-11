function [t,x,y,theta,cptim]=cablenl   
% [t,x,y,theta,cptim]=cablenl
% Example: cablenl
% ~~~~~~~~~~~~~~~~
% Numerical integration of the matrix 
% differential equations for the nonlinear 
% dynamics of a cable of rigid links with
% the outer ends of the cable fixed. 
%
% t     - time vector for the solution
% x,y   - matrices with nodal coordinates
%         stored in the columns. The time
%         history of point j is in x(:,j)
%         and y(:,j) 
% theta - matrix with inclination angles
%         stored in the columns
% cptim - number of seconds to integrate
%         the equations of motion
%
% User m functions required: 
%         plotmotn, equamo

clear all; close;

% Make variables global for use by 
% function equamo
global first_ n_ m_ len_ grav_ b_ mas_ py_

fprintf('\nNONLINEAR DYNAMICS OF A ')
fprintf('FALLING CABLE\n')
fprintf(...
'\nNote: The calculations take awhile\n')

% Set up data for a cable of n_ links, 
% initially arranged in a triangular 
% deflection configuration.

% parameter controlling initialization of 
% auxiliary parameters used in function 
% equamo
first_=1; 
% number of links in the cable
n_=16; n=n_; nh=n_/2;           
% vector of lengths and gravity constant
len_=1/n*ones(n,1); grav_=1;   
% vector of mass constants
m_=ones(1,n_)/n_;              

% initial position angles
th0=pi/4*[ones(nh,1);-ones(nh,1)]; 
td0=zeros(size(th0)); z0=[th0;td0];

% time limits, integration tolerances, 
% and the number of solution points
tmin=0; tmax=25; nt=201;
t=linspace(0,tmax,nt)';
tolrel=1e-6; tolabs=1e-8; len=len_;

% Perform the numerical integration using a 
% variable stepsize Runge-Kutta integrator
tic; 
odetol=odeset('RelTol',tolrel,'AbsTol',tolabs);
[t,w]=ode45(@equamo,t,z0,odetol);
theta=w(:,1:n); cptim=toc; 

% Compute node point coordinates
Z=[zeros(nt,1),repmat(len',nt,1).*exp(i*theta)];
Z=cumsum(Z.').'; x=real(Z); y=imag(Z); 

% Plot the horizontal position of the midpoint
clf; plot(t,x(:,1+n_ /2));  
ylabel('x coordinate'); xlabel('time')
title(['Horizontal Position of the ' ...
       'Cable Midpoint'])
grid on; figure(gcf); 
% print -deps xmidl

disp(' '), disp(...
'Press [Enter] to see the error growth curve');
pause, close

% Show error growth indicated by symmetry 
% loss of the vertical deflection symmetry.
% An approximately linear trend on the semilog
% plot indicates exponential growth of the error.
unsymer=sqrt(sum((y-y(:,end:-1:1)).^2,2));
hold off; axis('normal'); clf;
semilogy(t,unsymer); 
xlabel('time'); ylabel('asymmetry error');
title(['Growing Loss of Symmetry in ' ...
       'Vertical Deflection']);
grid on; figure(gcf);
% print -deps unsymerr

disp(' '), disp(...
'Press [Enter] to see the response animation');

% Show animation of the cable response 
disp(' ')
disp('The  motion can be animated or a trace')
disp('can be shown for successive positions')
disp(['between t = ',num2str(tmin),...
      ' and t = ',num2str(tmax)])

% Plot the position for different times limits
titl='CABLE MOTION FOR T = ';
while 1
  disp(' '), disp(...
  ['Choose a plot option (1 <=> animate, ',...
   ' 2 <=> trace,'])
  opt=input('3 <=> stop)  > ? ');
  if opt==3, break, end
  disp(' '), disp(...
  'Give a time vector such as 0:.1:15')
  Tp=input('Time vector > ? ','s'); 
  if isempty(Tp), break, end
  tp=eval(Tp); tp=tp(:); T=[titl,Tp];
  xp=interp1q(t,x,tp); yp=interp1q(t,y,tp);
  if opt ==1, plotmotn(xp,yp,T)
  else, plotmotn(xp,yp,T,1), end
end
fprintf('\nAll Done\n')

%=============================================

function plotmotn(x,y,titl,isave)
%
% plotmotn(x,y,titl,isave)
% ~~~~~~~~~~~~~~~~~~~~
% This function plots the cable time 
% history described by coordinate values 
% stored in the rows of matrices x and y.
%
% x,y   - matrices having successive rows 
%         which describe position 
%         configurations for the cable
% titl  - a title shown on the plots
% isave - parameter controlling the form 
%         of output. When isave is not input, 
%         only one position at a time is shown
%         in rapid succession to animate the
%         motion. If isave is given a value,
%         then successive are all shown at
%         once to illustrate a kinematic 
%         trace of the motion history.
%
% User m functions called:  none
%----------------------------------------------

% Set a square window to contain all 
% possible positions
[nt,n]=size(x); 
if nargin==4, save =1; else, save=0; end
xmin=min(x(:)); xmax=max(x(:));
ymin=min(y(:)); ymax=max(y(:)); 
w=max(xmax-xmin,ymax-ymin)/2;
xmd=(xmin+xmax)/2; ymd=(ymin+ymax)/2;  
hold off; clf; axis('normal'); axis('equal'); 
range=[xmd-w,xmd+w,ymd-w,ymd+w];
title(titl)
xlabel('x axis'); ylabel('y axis')
if save==0
  for j=1:nt
    xj=x(j,:); yj=y(j,:);
    plot(xj,yj,'-k',xj,yj,'ok');
    axis(range), axis off
    title(titl)
    figure(gcf), drawnow, pause(.1)
  end
  pause(2)
else
  hold off; close
  for j=1:nt
    xj=x(j,:); yj=y(j,:);
    plot(xj,yj,'-k',xj,yj,'ok');
    axis(range), axis off, hold on
  end
  title(titl)
  figure(gcf), drawnow, hold off, pause(2)
end

% Save plot history for subsequent printing 
% print -deps plotmotn

%=============================================

function zdot=equamo(t,z)
%
% zdot=equamo(t,z)
% ~~~~~~~~~~~~~~~~
% Equation of motion for a cable fixed at 
% both ends and loaded by gravity forces only
%
% t        current time value
% z        column vector defined by
%          [thet(t);theta'(t)]
% zdot     column vector defined by 
%          the concatenation
%          z'(t) = [theta'(t);theta''(t)]
%
% User m functions called:  none.
%----------------------------------------------

% Values accessed as global variables
global first_ n_ m_ len_ grav_ b_ mas_ py_

% Initialize parameters first time 
% function is called
if first_==1, first_=0;
% mass parameters
  m_=m_(:); b_=flipud(cumsum(flipud(m_))); 
  mas_=b_(:,ones(n_,1)); 
  mas_=tril(mas_)+tril(mas_,-1)';
% load effects from gravity forces
  py_=-grav_*(b_-b_(n_)); 
end

% Solve for zdot = [theta'(t); theta''(t)];
n=n_; len=len_; 
th=z(1:n); td=z(n+1:2*n); td2=td.*td;
x=len.*cos(th); y=len.*sin(th);

% Matrix of mass coefficients and 
% constraint conditions
amat=[[mas_.*(x*x'+y*y'),x,y];
      [x,y;zeros(2,2)]'];

% Right side vector involves applied forces 
% and inertial terms
bmat=x*y'; bmat=mas_.*(bmat-bmat');

% Solve for angular acceleration. 
% Most computation occurs here.
soln=amat\[x.*py_+bmat*td2; y'*td2; -x'*td2];

% Final result needed for use by the  
% numerical integrator
zdot=[td; soln(1:n)];
