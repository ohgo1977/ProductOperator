% Exp090_refocusedINEPT_PhaseCycling.m
% refocused INEPT I => S
% Example to check phase cycling.
% Keeler, J., Understanding NMR Spectroscopy (1st Ed.), Wiley, 2005.
% pp. 174 - 175.

clear
close all

% % 2-steps
% phid = 1:2;
% ph1tab = [0,2];% I 90
% ph2tab = [0];  % S INEPT 1st 180
% ph3tab = [0];  % I INEPT 1st 180
% ph4tab = [0];  % S INEPT 2nd 90
% ph5tab = [1];  % I INEPT 2nd 90
% ph6tab = [0];  % S INEPT 3rd 180
% ph7tab = [0];  % I INEPT 3rd 180
% phRtab = [0,2];% Receiver

% 16-steps
phid = 1:16;
ph1tab = [0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2];% I 90
ph2tab = [0,2,0,2];                        % S INEPT 1st 180
ph3tab = [0,2,0,2];                        % I INEPT 1st 180
ph4tab = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3];% S INEPT 2nd 90
ph5tab = [1,1,3,3];                        % I INEPT 2nd 90
ph6tab = [0,2,0,2,1,3,1,3];                % S INEPT 3rd 180
ph7tab = [0,2,0,2];                        % I INEPT 3rd 180
phRtab = [0,0,2,2,1,1,3,3];                % Receiver

%
syms B J t1 t2

% Initial State
rho_ini = PO(2,{'Iz' 'Sz'},{B 1});


% IS system
a0_M = [];
rho_M = [];
rho_total = 0;
for ii = phid
    fprintf(1,'\nii: %2d\n',ii);            
    ph1 = PO.phmod(ph1tab,ii);
    ph2 = PO.phmod(ph2tab,ii);
    ph3 = PO.phmod(ph3tab,ii);
    ph4 = PO.phmod(ph4tab,ii);
    ph5 = PO.phmod(ph5tab,ii);
    ph6 = PO.phmod(ph6tab,ii);
    ph7 = PO.phmod(ph7tab,ii);
    phR = PO.phmod(phRtab,ii);

    % OOP dot-style, CS ommitted, Pulse positions moved.
    rho = rho_ini;                                          % Preparation of the initial rho  
    rho.dispPOtxt();
    rho = rho.pulse('I',ph1,1/2*pi);                        % I 90 pulse
    rho = rho.simpulse({'I' 'S'},{ph3 ph2},{pi pi});        % I,S 180 pulses
    rho = rho.jc('IS',pi*J*2*t1);                           % J-coupling evolution
    rho = rho.simpulse({'I' 'S'},{ph5 ph4},{1/2*pi 1/2*pi});% I,S 90 pulses
    rho = rho.simpulse({'I' 'S'},{ph7 ph6},{pi pi});        % I,S 180 pulses
    rho = rho.jc('IS',pi*J*2*t2);                           % J-coupling evolution

    rho_detect = receiver(rho,phR);
    rho_total = rho_detect + rho_total;

    [a0_V, rho_V] = rho.SigAmp({'S'},phR); % Detection
    a0_M = cat(1,a0_M,a0_V);
    rho_M = cat(1,rho_M,rho_V);
end
rho_final = observable(rho_total,{'S'});