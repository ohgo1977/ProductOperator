clear
close all

spin_label_cell = {'I1' 'I2'};
PO.create(spin_label_cell);
gH_cell = PO.v2cell(gH,spin_label_cell);

for ii = 1:1
    % if ii ==1
    %     rho = (exp(-t1*(o1 + J12*pi)*1i)*1i)/2*I1p*I2a;
    % elseif ii == 2
    %     rho = (exp(-t1*(o1 - J12*pi)*1i)*1i)/2*I1p*I2b;
    % elseif ii == 3
    %     rho = -(exp(t1*(o1 + J12*pi)*1i)*1i)/2*I1m*I2a;
    % elseif ii == 4
    %     rho = -(exp(t1*(o1 - J12*pi)*1i)*1i)/2*I1m*I2b;
    % end
    if ii ==1
        rho = I1p*I2a;
    elseif ii == 2
        rho = I1p*I2b;
    elseif ii == 3
        rho = I1m*I2a;
    elseif ii == 4
        rho = I1m*I2b;
    end
    rho.disp = 1;
    rho = pfg(rho, G, gH_cell);
    rho = rho.pulse({'I*'},{'x'},{1/2*pi});
    rho = pfg(rho, G, gH_cell);
    dispPO(rho);
    rho_dephase = dephase(rho);
    dispPO(rho_dephase);
    if ii == 1
        rho_total = rho;
        rho_dephase_total = rho_dephase;
    else
        rho_total = rho + rho_total;
        rho_dephase_total = rho_dephase + rho_dephase_total;
    end
end