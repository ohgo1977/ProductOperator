% Exp150_RefocusingPulse_PFG.m
% Keeler, J., Understanding NMR Spectroscopy, p. 406, 11.12.3
% Gradient G - 180+d pulse - Gradient G 
% The selection of p => -p pathway.
% "Cleaning up" the results of an imperfect 180 pulse.

clear
close all

syms G gH d
pfg_switch = 1;

ini_status = 'DQ';
switch ini_status
    case 'SQ'
        spin_label_cell = {'I1'};
        rho = PO(1,{'I1p'},{1},spin_label_cell);% SQ
    case 'DQ'
        spin_label_cell = {'I1' 'I2'};
        rho = PO(2,{'I1pI2p'},{1},spin_label_cell);% DQ
    case 'TQ'
        spin_label_cell = {'I1' 'I2' 'I3'};
        rho = PO(3,{'I1pI2pI3p'},{1},spin_label_cell);% TQ
end
% % Alternative way to create rho from spin_label_cell
% ns = length(spin_label_cell);
% M_in = zeros(2^ns,2^ns);
% M_in(1,end) = 1;% I1pI2p...Inp
% rho = PO.M2pol(M_in,spin_label_cell);% Speed is a bit slower than PO().

dispPOtxt(rho);
gH_cell = PO.v2cell(gH,spin_label_cell);

% PFG
if pfg_switch == 1
    rho = pfg(rho, G, gH_cell);
end

% Imperfect 180 pulse (pi + d)
rho = simpulse(rho,{'I*'},{'x'},{pi + d});

% PFG
if pfg_switch == 1
    rho = pfg(rho, G, gH_cell);
end

dispPO(rho);

rho = dephase(rho);
dispPO(rho);