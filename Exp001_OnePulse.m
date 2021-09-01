clear
close all
% Examples of input parameters for PO.pulse().

rho = PO(1,{'Iz'});% Initial State
rho.dispPOtxt();

% All four cases below are equivalent.
pulse_writing = 'exp4';
switch pulse_writing
    case 'exp1'
        rho = rho.pulse('I','x',1/2*pi);% I90x-pulse
    case 'exp2'
        rho = rho.pulse(1,'x',1/2*pi);  % I90x-pulse
    case 'exp3'
        rho = rho.pulse('I',0,1/2*pi);  % I90x-pulse
    case 'exp4'
        rho = rho.pulse(1,0,1/2*pi);    % I90x-pulse
end