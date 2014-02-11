% Program svingetid
% Beregner svingetiden for en enkel, matematisk pendel
% fra den matematiske løsningen for store utslag
clear
thet0g = (5 : 5: 175);
n = length(thet0g);
thet0 = thet0g*pi/180;% theta(0) i radianer
T = zeros(n,1);
% T = fjerdeparten av perioden ( dimensjonsløs tid)
% for store utslag.
for l = 1:n
    k = sin(thet0(l)/2);
    k2 = k*k;
    T(l) = ellipke(k2);
end
% Plotting av resultater
plot(thet0g,T)
grid on
FS = 'FontSize';
% st = sprintf('\\theta_0 = %4.2f\\circ',thet0g);
title('Svingetiden for en kvart-period',FS,13);
xlabel('\theta_{0}(grader)',FS,14);
ylabel('T','Rotation',0,FS,13);
shg