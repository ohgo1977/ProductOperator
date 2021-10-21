% Exp090_refocusedINEPT_PhaseCycling_PS.m

% Para begin %
phid = 1:2;
ph_cell{1} = [0,2];% ph1
ph_cell{2} = [0];  % ph2
ph_cell{3} = [0];  % ph3
ph_cell{4} = [0];  % ph4
ph_cell{5} = [1];  % ph5
ph_cell{6} = [0];  % ph6
ph_cell{7} = [0];  % ph7
phRtab = [0,2];    % phR
spin_label_cell = {'I' 'S'};
coef_cell = {}; % Special sym coefs
rho_ini = Iz*B + Sz;
obs_cell = {'S'};
% Para end %

% PS begin %
rho = rho.pulse('I',ph1,1/2*pi);                        
rho = rho.simpulse({'I' 'S'},{ph3 ph2},{pi pi});        
rho = rho.jc('IS',pi*J*2*t1);                           
rho = rho.simpulse({'I' 'S'},{ph5 ph4},{1/2*pi 1/2*pi});
rho = rho.simpulse({'I' 'S'},{ph7 ph6},{pi pi});        
rho = rho.jc('IS',pi*J*2*t2);                           
% PS end %