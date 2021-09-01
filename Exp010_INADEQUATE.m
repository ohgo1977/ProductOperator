clear
close all
% Levitt, M. H., Spin Dynamics(2nd Ed.), p.433.
% 2D-INADEQUATE using -45 deg phase shift

syms oI oS t
% OOP dot-style
rho = PO(2,{'Iz' 'Sz'});
rho.dispPOtxt();

States = 'cos';

switch States
    case 'cos'
        phi = 0;

   case 'sin'
        phi = -1/4*pi;
        
end

rho = rho.simpulse_phshift({'I' 'S'},{phi phi},{3/2*pi 3/2*pi});
rho = rho.jc('IS',pi/2);
rho = rho.simpulse_phshift({'I' 'S'},{phi phi},{1/2*pi 1/2*pi});     
rho = rho.simcs({'I' 'S'},{oI*t oS*t});
rho = rho.simpulse({'I' 'S'},{'y' 'y'},{pi/2 pi/2});

phR = 0;
[a0_V,rho_V] = rho.SigAmp('IS',phR);