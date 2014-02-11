% program btabell
% Kragbjelke med stor utbøyning.
% Beregner lastfaktor alfa, nedbøyning delta
% samt horisontal lengde lh som funksjon av
% theta0 = t0g i grader
clear
t0g = 56.4946093;
fprintf('theta0   k       u1      alfa     delta      lh \n')
    t0 = t0g*pi/180;
    st0 = sin(t0);
    k = sqrt(0.5*(1 + sin(t0)));
    u1 = asin(1/(sqrt(2)*k));
    [K,Fu1] = ellipfk(u1,k);
    alfa = K - Fu1;
    [E,Eu1] = ellipek(u1,k);
    d = 1 - 2*(E - Eu1)/alfa;
    lh = sqrt(2*sin(t0))/alfa;
    fprintf(' %5.2f %8.5f %8.5f %8.5f %8.5f %8.5f \n',t0g,k,u1,alfa,d,lh);