clear
close all
% Example of chemical shift evolution

% Symbolic constant
syms q

rho = PO(1,{'Iz'});% Initial State
rho.dispPOtxt();
rho = rho.pulse('I','y',1/2*pi);% I90y-pulse
rho = rho.cs('I',q);% CS evolution
