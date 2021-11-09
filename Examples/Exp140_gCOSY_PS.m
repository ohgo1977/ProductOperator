% Exp140_gCOSY_PS.m
% Berger, S.; Braun, S. 200 and More NMR Experiments A Practical Course, p. 526

% Para begin %
ph_cell{1} = [0 2];
ph_cell{2} = [0 0 2 2];
phRtab = [0 2];% Receiver
spin_label_cell = {'I1' 'I2'};
rho_ini = I1z; 
obs_cell = {'I*'};
phid = 1:4;
coef_cell = {}; % Special sym coefs
disp_bin = 1;
syms gH;
gH_cell = PO.v2cell(gH,spin_label_cell);
% Para end %

% PS begin %
rho = rho.pulse({'I*'},{ph1},{1/2*pi});
rho = rho.cs({'I1' 'I2'},{o1*t1 o2*t1});
rho = rho.jc({'I1I2'},{pi*J12*t1});
rho = pfg(rho, G, gH_cell);
rho = rho.pulse({'I*'},{ph2},{1/2*pi});
rho = pfg(rho, G, gH_cell);
rho = dephase(rho);
% PS end %