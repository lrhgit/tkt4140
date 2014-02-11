% program sorfig
% tegner fig.7.11
clear
omega = (1: 0.05: 1.95)';
iter = zeros(20,4);
iter(1,1) = 24; iter(1,2) = 83; iter(1,3) = 275; iter(1,4) = 880;
iter(2,1) = 22; iter(2,2) = 76; iter(2,3) = 252; iter(2,4) = 810;
iter(3,1) = 20; iter(3,2) = 69; iter(3,3) = 231; iter(3,4) = 745;
iter(4,1) = 17; iter(4,2) = 63; iter(4,3) = 212; iter(4,4) = 685;
iter(5,1) = 15; iter(5,2) = 58; iter(5,3) = 194; iter(5,4) = 629;
iter(6,1) = 11; iter(6,2) = 52; iter(6,3) = 177; iter(6,4) = 576;
iter(7,1) = 12; iter(7,2) = 47; iter(7,3) = 161; iter(7,4) = 526;
iter(8,1) = 13; iter(8,2) = 42; iter(8,3) = 146; iter(8,4) = 479;
iter(9,1) = 15; iter(9,2) = 37; iter(9,3) = 132; iter(9,4) = 434;
iter(10,1) = 17; iter(10,2) = 31; iter(10,3) = 118; iter(10,4) = 392;
iter(11,1) = 20; iter(11,2) = 25; iter(11,3) = 105; iter(11,4) = 351;
iter(12,1) = 22; iter(12,2) = 24; iter(12,3) = 92;  iter(12,4) = 312;
iter(13,1) = 25; iter(13,2) = 27; iter(13,3) = 79;  iter(13,4) = 274;
iter(14,1) = 29; iter(14,2) = 31; iter(14,3) = 66;  iter(14,4) = 238;
iter(15,1) = 35; iter(15,2) = 37; iter(15,3) = 52;  iter(15,4) = 202;
iter(16,1) = 45; iter(16,2) = 45; iter(16,3) = 47;  iter(16,4) = 166;
iter(17,1) = 55; iter(17,2) = 56; iter(17,3) = 59;  iter(17,4) = 129;
iter(18,1) = 75; iter(18,2) = 75; iter(18,3) = 80;  iter(18,4) =  85;
iter(19,1) = 113; iter(19,2) = 113; iter(19,3) = 114;  iter(19,4) = 118;
iter(20,1) = 233; iter(20,2) = 228; iter(20,3) = 226;  iter(20,4) = 221;
plot(omega,iter,'k')
grid on
ylim([0 300])
xlabel('\omega','FontSize',14)
ylabel('iterasjoner','FontSize',14)
shg

