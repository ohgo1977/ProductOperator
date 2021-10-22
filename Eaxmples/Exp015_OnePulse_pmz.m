% Exp015_OnePulse_pmz.m
% One pulse experiment in the lowering/raising operator basis.

clear
close all


rho = PO(2,{'Iz'});% Initial State
rho = xyz2pmz(rho);% Since the default basis of Iz is xyz, it is necessary to convert the basis to pmz.
rho.dispPOtxt();
rho = rho.pulse('I','y',1/2*pi);% I90x-pulse
