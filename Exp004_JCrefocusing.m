clear
close all
% Keeler, J., Understanding NMR Spectroscopy (1st Ed.), Wiley, 2005.
% pp. 168, Fig. 7.14
% I: t/2-   -t/2 => cs is not refocused
% S: t/2-180-t/2 => cs is refocused
%                   jc is refocused

syms J12 t oI oS
rho = PO(2,{'Ix' 'Sx'});% Initial State
rho.dispPOtxt();
rho = rho.cs('I',oI*t/2);
rho = rho.cs('S',oS*t/2);
rho = rho.jc('IS',pi*J12*t/2);

rho = rho.pulse('S','x',pi);% Refocusing pulse on S
% rho = rho.pulse('I','x',pi);% What if refocusing pulse is also applied to I.

rho = rho.cs('I',oI*t/2);
rho = rho.cs('S',oS*t/2);
rho = rho.jc('IS',pi*J12*t/2);
