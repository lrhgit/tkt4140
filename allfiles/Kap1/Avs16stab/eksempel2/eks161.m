% PROGRAM Eks161
% L�ser ligning i avsnitt 1.6.1 i kompendiet
% (Kompendium �rgang 2006 og tidligere)
% Eksempel p� et stivt systems
% 
clear
options = odeset('AbsTol', 1.0e-10,'Jacobian','on','JConstant','on',...
   'Events','on');
tspan = [0:0.005:0.2];
   y0 = [0.0 0.0 ];
   [t,y,te,ye,ie] = ode15s('fcn161',tspan,y0,options);
   fprintf('\   t      y       v\n');
   fprintf('  %12.4e  %12.5e  %12.5e\n',[t y]');
   %plot(t,y(:,1))
   %plot(te,ye(:,2))
   
   
 