% program btabell
% Kragbjelke med stor utbøyning.
% Beregner lastfaktor alfa, nedbøyning delta
% samt horisontal lengde lh som funksjon av
% theta0 = t0g i grader
clear
t0g = (0:5:80)';
t0g = [t0g;80.5;81;81.5;82;82.5;83;83.5;84;84.5;85;85.5;86.0;86.5;87;87.5;88;88.5;89];
n = length(t0g);
fprintf('theta0   k       u1      alfa     delta      lh \n')
for l = 1:n
    t0 = t0g(l)*pi/180;
    st0 = sin(t0);
    k = sqrt(0.5*(1 + sin(t0)));
    u1 = asin(1/(sqrt(2)*k));
    [K,Fu1] = ellipfk(u1,k);
    alfa = K - Fu1;
    [E,Eu1] = ellipek(u1,k);
    d = 1 - 2*(E - Eu1)/alfa;
    lh = sqrt(2*sin(t0))/alfa;
    fprintf(' %4.1f %8.5f %8.5f %8.5f %8.5f %8.5f \n',t0g(l),k,u1,alfa,d,lh);
end