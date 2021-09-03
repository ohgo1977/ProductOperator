clear
close all
% Intensity calculation of refocused INEPT in InS system (n = 1,2 or 3)
% Levitt, M. H., Spin Dynamics (2nd Ed.), pp. 440 - 442, pp.488 - 491.

InS = 'I2S';
syms B J t
switch InS
    case 'IS'
        % IS system
        
        % syms BI BS
        % rho = PO(2,{'IxSz' 'Sz' '1'},{BI/4 -BS/4 sym(1/4)},{'I' 'S'});
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
        rho = simpulse(rho,{'I1' 'I2' 'S3'},{'x' 'x' 'x'},{3/2*pi 3/2*pi pi});
        rho = jc(rho,'I1S3',1/2*pi);
        rho = jc(rho,'I2S3',1/2*pi);

        rho = simpulse(rho,{'I1' 'I2' 'S3'},{'y' 'y' 'y'},{1/2*pi 1/2*pi 1/2*pi});
        rho = simpulse(rho,{'I1' 'I2' 'S3'},{'x' 'x' 'x'},{pi pi pi});
        rho = jc(rho,'I1S3',pi*J*t);
        rho = jc(rho,'I2S3',pi*J*t);
        a0_V = rho.SigAmp('S3','y');
        
    case 'I3S'
        % I3S system
        rho = PO(4,{'I1z' 'I2z' 'I3z' 'S4z'},{B B B sym(1)},{'I1' 'I2' 'I3' 'S4'});
        rho = simpulse(rho,{'I1' 'I2' 'I3' 'S4'},{'x' 'x' 'x' 'x'},{3/2*pi 3/2*pi 3/2*pi pi});
        rho = jc(rho,'I1S4',1/2*pi);
        rho = jc(rho,'I2S4',1/2*pi);
        rho = jc(rho,'I3S4',1/2*pi);

        rho = simpulse(rho,{'I1' 'I2' 'I3' 'S4'},{'y' 'y' 'y' 'y'},{1/2*pi 1/2*pi 1/2*pi 1/2*pi});
        rho = simpulse(rho,{'I1' 'I2' 'I3' 'S4'},{'x' 'x' 'x' 'x'},{pi pi pi pi});
        rho = jc(rho,'I1S4',pi*J*t);
        rho = jc(rho,'I2S4',pi*J*t);        
        rho = jc(rho,'I3S4',pi*J*t);
        a0_V = rho.SigAmp('S4','y');
end