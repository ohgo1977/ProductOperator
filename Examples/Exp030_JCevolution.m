% Exp030_JCevolution.m
% Example of J-coupling evolution

clear
close all

PO.create({'I' 'S'});% Preparation of PO objects and symbolic constants
rho = Ix;% Initial State
rho.dispPOtxt();
rho = rho.jc('IS',pi*J12*t);% J-coupling eovlution
