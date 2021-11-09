% Exp110_3QF_COSY_PS.m
% Guntert, P. et al., J. Magn. Reson. Ser. A, 101, 103-105, 1993.
% Guntert, P. Int. J. Quant. Chem., 106, 344-350, 2006.

% Para begin %
ph_cell{1} = sym([0:5]*pi/3);% I 90
phRtab = [0 2];% Receiver
spin_label_cell = {'I1' 'I2' 'I3'};
rho_ini = I1z; 
% rho_ini = PO(length(spin_label_cell),{'I1z'},{1},spin_label_cell);
obs_cell = {'I*'};
phid = 1:6;
coef_cell = {}; % Special sym coefs
disp_bin = 1;
% Para end %

% PS begin %
rho = rho.pulse_phshift({'I*'},{ph1},{1/2*pi});
rho = rho.cs({'I*'},{o1*t1});
rho = rho.jc({'I1I2' 'I1I3'},{pi*J12*t1 pi*J13*t1});
rho = rho.pulse_phshift({'I*'},{ph1},{1/2*pi});
rho = rho.pulse({'I*'},{0},{1/2*pi});
% PS end %