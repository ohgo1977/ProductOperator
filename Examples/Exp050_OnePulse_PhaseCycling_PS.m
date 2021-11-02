% Exp050_OnePulse_PhaseCycling_PS.m
% Example of writing a pulse sequence with phase cycling

% Para begin %
phid = 1:4;
ph_cell{1} = [1,2,3,0];     % Phase for 90-pulse
phRtab = [0,1,2,3];         % Receiver phase
spin_label_cell = {'I'};
coef_cell = {};             % Special sym coefs
rho_ini = Iz;
obs_cell = {'I'};
% Para end %

% PS begin %
rho = rho.pulse({1},{ph1},{pi/2});% 90-pulse                    
% PS end %
