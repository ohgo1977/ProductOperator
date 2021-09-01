clear
close all
% Keeler, J., Understanding NMR Spectroscopy (1st Ed.), Wiley, 2005.
% pp. 168, Fig. 7.14
% I: t/2-   -t/2 => cs is not refocused
% S: t/2-180-t/2 => cs is refocused
%                   js is refocused
syms J12 t oI oS
rho = PO(2,{'Ix' 'Sx'});% Initial State
rho.dispPOtxt();
rho = rho.cs('I',oI*t/2);
rho = rho.cs('S',oS*t/2);
rho = rho.jc('IS',pi*J12*t/2);

rho = rho.pulse('S','x',pi);% Refocusing pulse on S

rho = rho.cs('I',oI*t/2);
rho = rho.cs('S',oS*t/2);
rho = rho.jc('IS',pi*J12*t/2);
