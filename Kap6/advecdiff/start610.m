%program start610
% Tegner startprofilet (fig. 6.20a) som er en
% "flosshatt" (Top hat) for adveksjon-diffusjons-
% ligningen i avsnitt 6.10
clear; close
np = 100; % Antall punkt
np1 = np + 1;
dx = 1/np;
x = (0:dx:1)';
% Allokering
ui = zeros(np1,1); 
%
jc = round(np/2) + 1; % Midtpunkt
nc = jc - 1 ;
hw = round(np*0.1);
jstart = jc - hw;
jend = jc + hw;
% Startverdier 
for j = jstart:jend
    ui(j) = 0.5;
end
plot(x,ui);
grid
ylim([-0.1 0.6]);
FS = 'FontSize'; FW = 'FontWeight';
% % Overskrift for fig. 6.20a
title('Startprofil (flosshatt)',FS,14,FW,'Bold')

