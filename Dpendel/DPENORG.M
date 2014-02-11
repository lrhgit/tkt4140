% Program dpendel
% Løser bevegelsesligningene for en 
% dobbel-pendel med store utslag.
% Se Irgens : "Dynamikk", Kap. 9, eks.9.6
%
global g m1 m2 l1 l2;
% ==== Data ===
g = 9.81;
m1 = 0.5 ; m2 = 2.0;
l1 = 1.5; l2 = 1.0;
% === Variable ===
% Setter theta1 = y1 og theta2 = y3
% slik at d(theta1)/dt = d(y1)/dt = y2
% og d(theta2)/dt = d(y3)/dt = y4
%
tstart = 0.0; tstopp = 5.0;
antall = 10;
tspan = linspace(tstart,tstopp,antall);
% === Startverdier ====
theta1g = 90; % theta1 i grader
theta1 = theta1g*pi/180;
theta2g = 180;  % theta2 i grader
theta2 = theta2g*pi/180;
dtheta1 = 0.0; % Vinkelhastighet for l1
dtheta2 = 0.0; % Vinkelhastighet for l2
y0 =[theta1; dtheta1; theta2; dtheta2];
options = odeset('RelTol',1.0e-4,'Mass',@mass);
[t,y] = ode45(@fcndpend,tspan,y0,options);
x1 = l1*sin(y(:,1));
y1 = l1*cos(y(:,1));
x2 = x1 + l2*sin(y(:,3));
y2 = y1 + l2*cos(y(:,3));
x0 = 0.0; y0 = 0.0 ;  p1 = l1 + l2;
n = length(x1);
clf reset;
%hold on
for i = 1:n
    xl1 = [x0;x1(i)]; yl1 =[y0;y1(i)];
    xl2 = [x1(i);x2(i)]; yl2 = [y1(i); y2(i)];
    plot(xl1,yl1,xl2,yl2,'-ro','LineWidth',3);
    axis([-p1 p1 -p1 pl],'ij');
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    pause(0.1);
end
%hold off
% === Plotting ===
% figure('backingstore','off');
% axes('DrawMode','normal','box','on');
% axis([min(x1) max(x1) min(y1) max(y1)]);
% line_handle = line(x(1:2),y(1:2));
% set(line_handle,'LineWidth',2,'EraseMode','none');
% for counter = 1: length(x1)-1
%     set(line_handle,...
%         'XData',[x(counter + [0 +1])],...
%         'YData',[y(counter + [0 +1])]);
%     pause;
%     drawnow;
% end
%axis('ij');