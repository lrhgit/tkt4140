%=============== Program hekule =========================
% The program computes the velocity of a sphere
% falling vertically in a fluid. The equation 
% of motion is integrated using the Heun method.
%The first part output a table from t = 0 to t = tend.
% The second part integrates until the terminal velocity
% is attained and plots v versus t.
%=======================================================

nu = 1.5e-5 ;  % Kinematical viscosity [m^2/s]
rof = 1.22  ;  % Density of fluid [kg/m^3]
rol = 7850.0;  % Density of sphere [kg/m^3]
d   =  0.01 ;  % Diameter of sphere [m]
dt  = 0.1   ;  % Timestep [s]
tend = 2.0  ;  % Max. time [s]
fprintf(' Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
fprintf(' Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
fprintf(' Density of sphere ..... rol  = %10.3e kg/m^3 \n',rol);
fprintf(' Diameter of sphere .... d    = %10.3e m \n',d);
fprintf(' Timestep .............. dt   = %10.3e s \n',dt);
fprintf(' Max. time ............. tend = %10.3e s \n\n',tend);

g = 9.81    ;  % Gravity [N/kg]
ro = rof/rol;
A = 1.0 + 0.5*ro ;
B = (1.0 - ro)*g ;
C = 0.75*ro/d;
nsteps = round(tend/dt);
v = 0.0 ; t = 0.0;
fprintf('       t(s)       v(m/s)        Re \n\n');
% === OUTPUT A TABLE ===
for k = 1:nsteps
   t = k*dt;
   va = abs(v); Re = va*d/nu; CD = CDkule(Re);
   f = (B - C*v*va*CD)/A;
   vp = v + dt*f; % Predicted velocity
   vap = abs(vp); Re = vap*d/nu; CDp = CDkule(Re);
   fp = (B - C*vp*vap*CDp)/A;
   v = v + 0.5*dt*(f + fp); % Corrected velocity
   fprintf(' %10.2f  %10.3f %15.3e \n',t,v,Re);
end
% === COMPUTE UNTIL THE TERMINAL VELOCITY IS ATTAINED ===
% We collect v and t in vectors for plotting.
v(1) = 0.0 ; t(1) = 0.0; Re = 0.0; vt = 1.0;
epsi = 5.0e-3; test = 1; k = 0;
while test > epsi
    k = k + 1;
    t(k+1) = k*dt;
    va = abs(v(k)); Re = va*d/nu; CD = CDkule(Re);
    f = (B - C*v(k)*va*CD)/A;
    vp = v(k) + dt*f; % Predicted velocity
    vap = abs(vp); Re = vap*d/nu; CDp = CDkule(Re);
    fp = (B - C*vp*vap*CDp)/A;
    v(k+1) = v(k) + 0.5*dt*(f + fp); % Corrected velocity
    if k > 1
        vt = sqrt(B/(CD*C));
    else
        vt = 0;
    end
    test = abs((vt - v(k+1))/v(k+1));
end
fprintf('\n Terminal velocity = %7.3f (m/s) at t = %7.3f s \n',v(end),t(end));
plot(t,v);
FS = 'FontSize'; FW = 'FontWeight';
xlabel('t(s)',FS,14,FW,'Bold')
ylabel('v(m/s)',FS,14,FW,'Bold')
   


