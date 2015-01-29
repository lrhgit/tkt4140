% program sorfig3
% Generate fig.7.14
clear
Tstart = (0.5: 0.001: 0.52)';
% Computed optimal omega
omega = zeros(21,4);
omega(1,1)  = 1.14; omega(1,2)  = 1.41; omega(1,3)  = 1.64;  omega(1,4)  = 1.80;
omega(2,1)  = 1.38; omega(2,2)  = 1.59; omega(2,3)  = 1.74;  omega(2,4)  = 1.80;
omega(3,1)  = 1.42; omega(3,2)  = 1.62; omega(3,3)  = 1.77;  omega(3,4)  = 1.84;
omega(4,1)  = 1.44; omega(4,2)  = 1.64; omega(4,3)  = 1.78;  omega(4,4)  = 1.87;
omega(5,1)  = 1.44; omega(5,2)  = 1.65; omega(5,3)  = 1.79;  omega(5,4)  = 1.88;
omega(6,1)  = 1.45; omega(6,2)  = 1.66; omega(6,3)  = 1.79;  omega(6,4)  = 1.88;
omega(7,1)  = 1.46; omega(7,2)  = 1.66; omega(7,3)  = 1.80;  omega(7,4)  = 1.88;
omega(8,1)  = 1.47; omega(8,2)  = 1.66; omega(8,3)  = 1.80;  omega(8,4)  = 1.89;
omega(9,1)  = 1.48; omega(9,2)  = 1.66; omega(9,3)  = 1.80;  omega(9,4)  = 1.89;
omega(10,1) = 1.48; omega(10,2) = 1.67; omega(10,3) = 1.81;  omega(10,4) = 1.89;
omega(11,1) = 1.48; omega(11,2) = 1.67; omega(11,3) = 1.81;  omega(11,4) = 1.89;
omega(12,1) = 1.48; omega(12,2) = 1.67; omega(12,3) = 1.81;  omega(12,4) = 1.89;
omega(13,1) = 1.48; omega(13,2) = 1.67; omega(13,3) = 1.81;  omega(13,4) = 1.90;
omega(14,1) = 1.48; omega(14,2) = 1.67; omega(14,3) = 1.81;  omega(14,4) = 1.90;
omega(15,1) = 1.48; omega(15,2) = 1.67; omega(15,3) = 1.81;  omega(15,4) = 1.90;
omega(16,1) = 1.48; omega(16,2) = 1.68; omega(16,3) = 1.81;  omega(16,4) = 1.90;
omega(17,1) = 1.48; omega(17,2) = 1.68; omega(17,3) = 1.81;  omega(17,4) = 1.90;
omega(18,1) = 1.48; omega(18,2) = 1.68; omega(18,3) = 1.81;  omega(18,4) = 1.90;
omega(19,1) = 1.48; omega(19,2) = 1.68; omega(19,3) = 1.81;  omega(19,4) = 1.90;
omega(20,1) = 1.48; omega(20,2) = 1.68; omega(20,3) = 1.81;  omega(20,4) = 1.90;
omega(21,1) = 1.49; omega(21,2) = 1.68; omega(21,3) = 1.81;  omega(21,4) = 1.90;

plot(Tstart,omega,'k')
grid on
xlim([0.5 0.52])
ylim([1.0 2.0])
xlabel('startverdier','FontSize',14)
ylabel('\omega_{opt}','FontSize',14,'FontWeight','Bold')
shg

