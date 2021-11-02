% Exp020_CSevolution.m
% Example of chemical shift evolution

clear
close all

% Symbolic constant
syms q

rho = PO(1,{'Iz'});% Initial State
rho.dispPOtxt();
rho = rho.pulse({'I'},{'y'},{1/2*pi});% I90y-pulse
rho = rho.cs({'I'},{q});% CS evolution
