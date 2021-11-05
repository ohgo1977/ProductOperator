% Exp130_PFG_Phase.m
% Keeler, J., Understanding NMR Spectroscopy, p. 399, 11.11.1
% Calculation of spatially dependent phase created by PFG for coherence order P.

clear
close all

syms G gH

ini_status = 'DQ';
switch ini_status
    case 'SQ'
        spin_label_cell = {'I1'};
    case 'DQ'
        spin_label_cell = {'I1' 'I2'};
    case 'TQ'
        spin_label_cell = {'I1' 'I2' 'I3'};
end

ns = length(spin_label_cell);
gH_cell = PO.v2cell(gH,spin_label_cell);

for ii = 1:2^ns
    for jj = 1:2^ns
        M_in = zeros(2^ns,2^ns);
        if ii ~= jj
            M_in(ii,jj) = 1;
            rho = PO.M2pol(M_in,spin_label_cell);% Speed is a bit slower than PO().
            rho.disp = 0;
            rho = pfg(rho, G, gH_cell);
            dispPO(rho);
        end
    end
end