clear
close all
% Spin-echo (Hahn-echo) experiment with phase cycling.

% Phase tables
phid = 1:16;
ph1tab = [1,2,3,0];                         % Phase for 90-pulse
ph2tab = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]; % Phase for 180-pulse
phRtab = [0,3,2,1,2,1,0,3];                 % Receiver phase

% Symbolic constants
syms t oI d

% Initialization
a0_M = [];
rho_M = [];

% Pulse sequence with phase cycling
for ii = phid
    fprintf(1,'\nii: %2d\n',ii);            
    ph1 = PO.phmod(ph1tab,ii);
    ph2 = PO.phmod(ph2tab,ii);
    phR = PO.phmod(phRtab,ii);

    rho = PO(1,{'Iz'});% Initial State
    rho.dispPOtxt();
    rho = rho.pulse(1,ph1,1/2*pi);% 90-pulse
    rho = rho.cs(1,oI*t);% Chemical shift evolution
    rho = rho.pulse(1,ph2,pi+d);% 180+d-pulse, where d indicates the miscalibration of 180-pulse
 %   rho = rho.pulse(1,ph2,pi);% 180-pulse
    rho = rho.cs(1,oI*t);% Chemical shift evolution

    [a0_V,rho_V] = rho.SigAmp('I',phR);% Detection
    a0_M = cat(1,a0_M,a0_V);
    rho_M = cat(1,rho_M,rho_V);
end