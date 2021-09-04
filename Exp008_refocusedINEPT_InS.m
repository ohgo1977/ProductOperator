clear
close all
% Intensity calculation of refocused INEPT in InS system (n = 1,2 or 3)
% Levitt, M. H., Spin Dynamics (2nd Ed.), pp. 440 - 442, pp.488 - 491.

InS = 'I3S';
syms B J t
switch InS
    case 'IS'
        % IS system
        rho = PO(2,{'Iz' 'Sz'},{B sym(1)},{'I' 'S'});
        rho = simpulse(rho,{'I' 'S'},{'x' 'x'},{3/2*pi pi});
        rho = jc(rho,'IS',1/2*pi);
        
        rho = simpulse(rho,{'I' 'S'},{'y' 'y'},{1/2*pi 1/2*pi});
        rho = simpulse(rho,{'I' 'S'},{'x' 'x'},{pi pi});
        rho = jc(rho,'IS',pi*J*t);
        a0_V = rho.SigAmp('S','y');

    case 'I2S'
        % I2S system
        rho = PO(3,{'I1z' 'I2z' 'S3z'},{B B sym(1)},{'I1' 'I2' 'S3'});
        rho = simpulse(rho,{'I*' 'S3'},{'x' 'x'},{3/2*pi pi});
        rho = simjc(rho,{'I1S3' 'I2S3'},{1/2*pi 1/2*pi});

        rho = simpulse(rho,{'I*' 'S3'},{'y' 'y'},{1/2*pi 1/2*pi});
        rho = simpulse(rho,{'I*' 'S3'},{'x' 'x'},{pi pi});
        rho = simjc(rho,{'I1S3' 'I2S3'},{pi*J*t pi*J*t});
        a0_V = rho.SigAmp('S3','y');
        
    case 'I3S'
        % I3S system
        rho = PO(4,{'I1z' 'I2z' 'I3z' 'S4z'},{B B B sym(1)},{'I1' 'I2' 'I3' 'S4'});
        rho = simpulse(rho,{'I*' 'S4'},{'x' 'x'},{3/2*pi pi});
        rho = simjc(rho,{'I1S4' 'I2S4' 'I3S4'},{1/2*pi 1/2*pi 1/2*pi});

        rho = simpulse(rho,{'I*' 'S4'},{'y' 'y'},{1/2*pi 1/2*pi});
        rho = simpulse(rho,{'I*' 'S4'},{'x' 'x'},{pi pi});
        rho = simjc(rho,{'I1S4' 'I2S4' 'I3S4'},{pi*J*t pi*J*t pi*J*t});
        a0_V = rho.SigAmp('S4','y');
end