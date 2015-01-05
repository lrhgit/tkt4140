% Program test
% Diskret cosinus transformasjon av m ekvidistante punkt
% (xj,yj), j = 0,1,.. m-1 der x(j) = pi(j + 1/2)/m
% og m er et partall.
%
%        --- Input ---
% y: vektor med  m y-ordinater y0,y1,y2,..,yn
%         (kolonne eller linjevektor)
%        --- Output ---
% a0 : 1. ledd i cosinusrekke
% a  : vektor med koeffisienter for cosinusleddene.
%      Antall ledd er m/2 
% a 
%
m = 8;
i = 0: m - 1;
x = pi*(i + 1/2)/m;
y = x;
[a0,a] = dct(y);
