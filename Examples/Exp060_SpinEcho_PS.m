% Exp060_SpinEcho_PS.m
% Spin-echo (Hahn-echo) experiment with phase cycling.
% Effect of the miscalibration of 180 pulse can be checked.

% Para begin %
phid = 1:16;
ph_cell{1} = [1,2,3,0];                         % Phase for 90-pulse
ph_cell{2} = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]; % Phase for 180-pulse
phRtab = [0,3,2,1,2,1,0,3];    % phR
spin_label_cell = {'I'};
coef_cell = {}; % Special sym coefs
rho_ini = Iz;
obs_cell = {'I'};
% Para end %

% PS begin %
rho = rho.pulse({1},{ph1},{1/2*pi});% 90-pulse
rho = rho.cs({1},{oI*t});% Chemical shift evolution
% rho = rho.pulse({1},{ph2},{pi});% 180-pulse
rho = rho.pulse({1},{ph2},{pi+d});% 180+d-pulse, where d indicates the miscalibration of 180-pulse
rho = rho.cs({1},{oI*t});% Chemical shift evolution                       
% PS end %
