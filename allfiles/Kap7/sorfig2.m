% program sorfig2
% Generate fig.7.13
clear
Tstart = (0.5: 0.001: 0.52)';
% iterations for theoretical optimal omega
iter = zeros(21,4);
iter(1,1)  = 9;  iter(1,2)  = 17; iter(1,3)  = 34;  iter(1,4)  = 67;
iter(2,1)  = 18; iter(2,2)  = 30; iter(2,3)  = 46;  iter(2,4)  = 68;
iter(3,1)  = 21; iter(3,2)  = 37; iter(3,3)  = 59;  iter(3,4)  = 90;
iter(4,1)  = 23; iter(4,2)  = 41; iter(4,3)  = 67;  iter(4,4)  = 105;
iter(5,1)  = 24; iter(5,2)  = 44; iter(5,3)  = 73;  iter(5,4)  = 116;
iter(6,1)  = 25; iter(6,2)  = 46; iter(6,3)  = 77;  iter(6,4)  = 125;
iter(7,1)  = 26; iter(7,2)  = 48; iter(7,3)  = 80;  iter(7,4)  = 132;
iter(8,1)  = 27; iter(8,2)  = 49; iter(8,3)  = 83;  iter(8,4)  = 138;
iter(9,1)  = 28; iter(9,2)  = 50; iter(9,3)  = 86;  iter(9,4)  = 143;
iter(10,1) = 28; iter(10,2) = 51; iter(10,3) = 88;  iter(10,4) = 147;
iter(11,1) = 29; iter(11,2) = 52; iter(11,3) = 90;  iter(11,4) = 151;
iter(12,1) = 29; iter(12,2) = 53; iter(12,3) = 92;  iter(12,4) = 155;
iter(13,1) = 29; iter(13,2) = 54; iter(13,3) = 94;  iter(13,4) = 158;
iter(14,1) = 30; iter(14,2) = 55; iter(14,3) = 95;  iter(14,4) = 161;
iter(15,1) = 30; iter(15,2) = 56; iter(15,3) = 97;  iter(15,4) = 164;
iter(16,1) = 30; iter(16,2) = 56; iter(16,3) = 98;  iter(16,4) = 167;
iter(17,1) = 31; iter(17,2) = 57; iter(17,3) = 99;  iter(17,4) = 169;
iter(18,1) = 31; iter(18,2) = 57; iter(18,3) = 101; iter(18,4) = 172;
iter(19,1) = 31; iter(19,2) = 58; iter(19,3) = 102; iter(19,4) = 174;
iter(20,1) = 31; iter(20,2) = 59; iter(20,3) = 103; iter(20,4) = 176;
iter(21,1) = 32; iter(21,2) = 59; iter(21,3) = 104; iter(21,4) = 178;
% iterations for computed otimal omega
iter2 = zeros(21,4);
iter2(1,1)  = 8;  iter2(1,2)  = 16; iter2(1,3)  = 32;  iter2(1,4)  = 64;
iter2(2,1)  = 13; iter2(2,2)  = 24; iter2(2,3)  = 42;  iter2(2,4)  = 67;
iter2(3,1)  = 14; iter2(3,2)  = 26; iter2(3,3)  = 47;  iter2(3,4)  = 85;
iter2(4,1)  = 15; iter2(4,2)  = 27; iter2(4,3)  = 50;  iter2(4,4)  = 87;
iter2(5,1)  = 15; iter2(5,2)  = 28; iter2(5,3)  = 51;  iter2(5,4)  = 91;
iter2(6,1)  = 15; iter2(6,2)  = 29; iter2(6,3)  = 52;  iter2(6,4)  = 94;
iter2(7,1)  = 15; iter2(7,2)  = 29; iter2(7,3)  = 53;  iter2(7,4)  = 98;
iter2(8,1)  = 15; iter2(8,2)  = 29; iter2(8,3)  = 54;  iter2(8,4)  = 99;
iter2(9,1)  = 15; iter2(9,2)  = 29; iter2(9,3)  = 54;  iter2(9,4)  = 100;
iter2(10,1) = 16; iter2(10,2) = 30; iter2(10,3) = 55;  iter2(10,4) = 100;
iter2(11,1) = 16; iter2(11,2) = 30; iter2(11,3) = 55;  iter2(11,4) = 102;
iter2(12,1) = 16; iter2(12,2) = 30; iter2(12,3) = 55;  iter2(12,4) = 103;
iter2(13,1) = 16; iter2(13,2) = 30; iter2(13,3) = 56;  iter2(13,4) = 104;
iter2(14,1) = 16; iter2(14,2) = 30; iter2(14,3) = 56;  iter2(14,4) = 105;
iter2(15,1) = 16; iter2(15,2) = 30; iter2(15,3) = 56;  iter2(15,4) = 105;
iter2(16,1) = 16; iter2(16,2) = 31; iter2(16,3) = 57;  iter2(16,4) = 105;
iter2(17,1) = 16; iter2(17,2) = 31; iter2(17,3) = 57;  iter2(17,4) = 105;
iter2(18,1) = 16; iter2(18,2) = 31; iter2(18,3) = 57;  iter2(18,4) = 105;
iter2(19,1) = 16; iter2(19,2) = 31; iter2(19,3) = 57;  iter2(19,4) = 106;
iter2(20,1) = 16; iter2(20,2) = 31; iter2(20,3) = 58;  iter2(20,4) = 106;
iter2(21,1) = 16; iter2(21,2) = 31; iter2(21,3) = 58;  iter2(21,4) = 107;
plot(Tstart,iter,'k',Tstart,iter2,'k-.')
grid on
xlim([0.5 0.52])
xlabel('startverdier','FontSize',14)
ylabel('iterasjoner','FontSize',14)
shg

