% Exp100_INADEQUATE_PS.m
% Levitt, M. H., Spin Dynamics(2nd Ed.), p.433.
% 2D-INADEQUATE using -45 deg phase shift

% Para begin %
phid = 1:1;
phRtab = [0];
% spin_label_cell = {'I1' 'I2'};
% rho_ini = I1z + I2z;
spin_label_cell = {'I' 'S'};
rho_ini = Iz + Sz;
coef_cell = {}; % Special sym coefs
obs_cell = {1 2};
phi_vec = [0 -1/4*pi];
States = 'sin';
phi_id = [contains(States,'cos') contains(States,'sin')];% switch syntax
phi = phi_vec(phi_id ~= 0);% switch syntax
% Para end %

% PS begin %
rho = rho.simpulse_phshift({1 2},{phi phi},{3/2*pi 3/2*pi});
rho = rho.jc([1 2],pi/2);
rho = rho.simpulse_phshift({1 2},{phi phi},{1/2*pi 1/2*pi});     
rho = rho.simcs({1 2},{o1*t o2*t});
rho = rho.simpulse({1 2},{'y' 'y'},{pi/2 pi/2});
% PS end %
