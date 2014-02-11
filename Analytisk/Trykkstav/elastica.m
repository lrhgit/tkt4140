% program elastica
% Trykkstav med stor utbøyning.
% Plotter elastica for gitte verdier av theta0g (grader)
% Bruker funksjonen fcnelast

% Tilsvarer tabell 2-4, avsnitt 2.7 i
% Timoshenko & Gere : Theory of Elastic Stability,
% 2. utgave McGraw-Hill, 1961
clear
t0g = (20:20:160)';
t0g = [t0g;176];
n = 40; % Deler staven i n deler
m = length(t0g);
hold on;
xlim([0 1]);
for l = 1:m
    theta0g = t0g(l);
    [x,y] = fcnelast(theta0g,n);
    plot(y,x);
end
hold off
title('Elastica som funksjon av \theta_{0}','FontSize',14)