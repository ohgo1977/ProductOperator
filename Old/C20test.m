clear
close all

syms B q
rho = PO(4,{'I1z' 'I2z' 'I3z' 'S4z'},{B B B sym(1)},{'I1' 'I2' 'I3' 'S4'});
rho.disp = 0;

rho1 = simpulse(rho,{'*'},{'x'},{pi/2});
dispPOtxt(rho1);

%rho2 = simpulse(rho,{'I*' 'S4'},{'x' 'x'},{pi/2 pi/2});
%dispPOtxt(rho2);

%rho3 = simpulse_phshift(rho,{'*'},{0},{pi/2});
%dispPOtxt(rho3);

%rho4 = simpulse_phshift(rho,{'I*' 'S4'},{0 0},{pi/2 pi/2});
%dispPOtxt(rho4);

rho5 = simpulse_phshift(rho,{'I*' 'S*'},{0 0},{pi/2 pi/2});
dispPOtxt(rho5);

%rho = rho1;
%rho1 = simcs(rho,{'*'},{q});
%dispPOtxt(rho1);
%rho2 = simcs(rho,{'I*' 'S4'},{q q});
%dispPOtxt(rho2);
%rho3 = simcs(rho,{'I*' 'S*'},{q q});
%dispPOtxt(rho3);