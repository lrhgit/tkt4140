% program etabell
% Elastica - Trykkstav med stor utbøyning.
% Beregner lastfaktor alfa, kraftforhold P/Pk, utbøyning delta
% samt vertikal lengde lv som funksjon av
% theta0 = t0g i grader
% Tilsvarer tabell 2-4, avsnitt 2.7 i
% Timoshenko & Gere : Theory of Elastic Stability,
% 2. utgave McGraw-Hill, 1961
clear
t0g = (0:10:170)';
t0g = [t0g;175;176;177;178;179];
n = length(t0g);
fprintf(' theta0     k       P/Pk      alfa     delta      lv \n')
for l = 1:n
    t0 = t0g(l)*pi/180;
    k = sin(t0/2);
    [K,E] = ellipke(k^2);
    alfa = K;
    ppk = (2*alfa/pi)^2;
    d = 2*k/alfa;    
    lv = 2*E/alfa - 1;
    fprintf(' %6.1f  %8.5f  %8.5f %8.5f %8.5f  %8.5f \n',t0g(l),k,ppk,alfa,d,lv);
end