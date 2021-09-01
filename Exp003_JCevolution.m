clear
close all
% Example of J-coupling evolution

% Symbolic consntants
syms J12 t

rho = PO(2,{'Ix'});% Initial State
rho.dispPOtxt();
rho = rho.jc('IS',pi*J12*t);% J-coupling eovlution
