% Exp050_OnePulse_PhaseCycling.m

clear
close all
% Example of writing a pulse sequence with phase cycling

% Phase tables
phid = 1:4;
ph1tab = [1,2,3,0];     % Phase for 90-pulse
phRtab = [0,1,2,3];     % Receiver phase

rho_ini = PO(1,{'Iz'});% Initial State

% Initialization
a0_M = [];
rho_M = [];
rho_total = 0;

% Pulse sequence with phase cycling
for ii = phid
    fprintf(1,'\nii: %2d\n',ii);            
    ph1 = PO.phmod(ph1tab,ii);
    phR = PO.phmod(phRtab,ii);

    rho = rho_ini;
    rho.dispPOtxt();% Display Initial state
    rho = rho.pulse(1,ph1,pi/2);% 90-pulse

    rho_detect = receiver(rho,phR);
    rho_total = rho_detect + rho_total;

    [a0_V,rho_V] = rho.SigAmp({'I'},phR);% Detection
    a0_M = cat(1,a0_M,a0_V);
    rho_M = cat(1,rho_M,rho_V);
end
rho_final = observable(rho_total,{'I'});