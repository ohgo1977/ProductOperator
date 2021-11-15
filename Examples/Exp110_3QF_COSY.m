% Exp110_3QF_COSY.m
% Guntert, P. et al., J. Magn. Reson. Ser. A, 101, 103-105, 1993.
% Guntert, P. Int. J. Quant. Chem., 106, 344-350, 2006.

clear
close all

phid = 1:6;
ph1tab = sym([0:5]*pi/3);% I 90
phRtab = [0 2];% Receiver

% Initial State
rho_ini = PO(3,{'I1z'},{1},{'I1' 'I2' 'I3'});
rho_ini.disp = 1;
PO.symcoef({'I1' 'I2' 'I3'})

a0_M = [];
rho_M = [];
rho_total = 0;

for ii = phid
    fprintf(1,'\nii: %2d\n',ii);            
    ph1 = PO.phmod(ph1tab,ii);
    phR = PO.phmod(phRtab,ii);

    rho = rho_ini;
    rho.dispPOtxt();

    rho = rho.pulse_phshift({'I*'},{ph1},{1/2*pi});
    rho = rho.cs({'I*'},{o1*t1});
    rho = rho.jc({'I1I2' 'I1I3'},{pi*J12*t1 pi*J13*t1});
    rho = rho.pulse_phshift({'I*'},{ph1},{1/2*pi});
    rho = rho.pulse({'I*'},{0},{1/2*pi});

    rho_detect = receiver(rho,phR);
    rho_total = rho_detect + rho_total;

    % [a0_V, rho_V] = rho.SigAmp({'I*'},phR); % Detection
    % a0_M = cat(1,a0_M,a0_V);
    % rho_M = cat(1,rho_M,rho_V);
end
rho_final = observable(rho_total,{'I*'});