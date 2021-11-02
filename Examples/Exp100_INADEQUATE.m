% Exp100_INADEQUATE.m
% Levitt, M. H., Spin Dynamics(2nd Ed.), p.433.
% 2D-INADEQUATE using -45 deg phase shift

clear
close all

syms oI oS t
rho = PO(2,{'Iz' 'Sz'});
% rho = xyz2pmz(rho);% Check the result in the pmz basis.
rho.dispPOtxt();

States = 'sin';

switch States
    case 'cos'
        phi = 0;

   case 'sin'
        phi = -1/4*pi;
        
end

rho = rho.pulse_phshift({'I' 'S'},{phi phi},{3/2*pi 3/2*pi});
rho = rho.jc({'IS'},{pi/2});
rho = rho.pulse_phshift({'I' 'S'},{phi phi},{1/2*pi 1/2*pi});     
rho = rho.cs({'I' 'S'},{oI*t oS*t});
rho = rho.pulse({'I' 'S'},{'y' 'y'},{pi/2 pi/2});

phR = 0;
[a0_V,rho_V] = rho.SigAmp({'I' 'S'},phR);