% Exp080_refocusedINEPT_InS.m

clear
close all
% Intensity calculation of refocused INEPT in InS system (n = 1,2 or 3)
% Levitt, M. H., Spin Dynamics (2nd Ed.), pp. 440 - 442, pp.488 - 491.

InS = 'I3S';
switch InS
    case 'IS'
        % IS system
        PO.create({'I1' 'S2'})
        rho = I1z*B + S2z;
        jc_cell = {'I1S2'};

    case 'I2S'
        % I2S system
        PO.create({'I1' 'I2' 'S3'});
        rho = I1z*B + I2z*B + S3z;
        jc_cell = {'I1S3' 'I2S3'};

    case 'I3S'
        % I3S system
        PO.create({'I1' 'I2' 'I3' 'S4'});
        rho = I1z*B + I2z*B + I3z*B + S4z;
        jc_cell = {'I1S4' 'I2S4' 'I3S4'};

end
q1 = 1/2*pi;
q1_cell = PO.v2cell(q1,jc_cell);

q2 = pi*J*t;
q2_cell = PO.v2cell(q2,jc_cell);

dispPOtxt(rho);
rho = simpulse(rho,{'I*' 'S*'},{'x' 'x'},{3/2*pi pi});
rho = simjc(rho,jc_cell,q1_cell);

rho = simpulse(rho,{'I*' 'S*'},{'y' 'y'},{1/2*pi 1/2*pi});
rho = simpulse(rho,{'I*' 'S*'},{'x' 'x'},{pi pi});
rho = simjc(rho,jc_cell,q2_cell);

rho_detect = receiver(rho,'x');
rho_final = observable(rho_detect,{'S*'});
dispPO(rho_final);
[a0_V,rho_V] = rho.SigAmp({'S*'},'x');