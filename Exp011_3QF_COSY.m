% clear
close all
% Guntert, P.

phid = 1:6;
ph1tab = sym([0:5]*pi/3);% I 90
phRtab = [0 2];% Receiver

rho_ini = PO(3,{'I1z'},{sym(1)},{'I1' 'I2' 'I3'});
rho_ini.disp = 1;

syms o1 t1 J12 J13

a0_M = [];
rho_M = [];
for ii = phid
    fprintf(1,'\nii: %2d\n',ii);            
    ph1 = PO.phmod(ph1tab,ii);
    phR = PO.phmod(phRtab,ii);

    rho = rho_ini;
    rho.dispPOtxt();

    rho = rho.simpulse_phshift({'I*'},{ph1},{1/2*pi});
    rho = rho.simcs({'I*'},{o1*t1});
    rho = rho.simjc({'I1I2' 'I1I3'},{pi*J12*t1 pi*J13*t1});
    rho = rho.simpulse_phshift({'I*'},{ph1},{1/2*pi});
    rho = rho.simpulse({'I*'},{0},{1/2*pi});

    dispPO(rho);
    [a0_V, rho_V] = rho.SigAmp2({'I*'},phR); % Detection

    % rho.disp = 1;
    % rho = rho.cs('I1',-phR*pi/2);
    % [a0_V, rho_V] = rho.SigAmp2({'I1'},0); % Detection

    a0_M = cat(1,a0_M,a0_V);
    rho_M = cat(1,rho_M,rho_V);
end