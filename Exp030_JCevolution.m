% Exp030_JCevolution.m

clear
close all
% Example of J-coupling evolution

PO.create({'I' 'S'});% Preparation of PO objects and symbolic constants

rho = Ix;% Initial State
rho.dispPOtxt();
rho = rho.jc('IS',pi*J12*t);% J-coupling eovlution
