clear
close all
% Keeler, J., Understanding NMR Spectroscopy, p. 406, 11.12.3
% Gradient G - 180+d pulse - Gradient G 
% The selection of P => -p pathway.
% "Cleaning up" the results of an imperfect 180 pulse.

syms G gH d
pfg_switch = 1;

% rho = PO(1,{'Ip'});% Ip
% Try the case from 2p to -2p
rho = PO(2,{'I1pI2p'},{sym(1)},{'I1' 'I2'});% I1pI2p

if pfg_switch == 1
    rho = pfg(rho, G, {gH gH});
end

rho = simpulse(rho,{'I*'},{'x'},{pi + d});% 180+d pulse

if pfg_switch == 1
    rho = pfg(rho, G, {gH gH});
end

dispPO(rho);

rho = dephase(rho);
dispPO(rho);

% Result w/ PFG is
% Iz
%   C:-sin(d)*(cos(G*gH)*1i + sin(G*gH))
%   T:1
% Ip
%   C:[(cos(G*gH)*1i + sin(G*gH))^2/2, -(cos(G*gH)*1i + sin(G*gH))^2/2]
%   T:[cos(d), 1]
% Im
%   C:[1/2, 1/2]
%   T:[cos(d), 1]
% Only Im is PFG-independent and Ip will dissaper by PFG.
%
% Result w/o PFG is
% Iz
%   C:-sin(d)*1i
%   T:1
% Ip
%   C:[-1/2, 1/2]
%   T:[cos(d), 1]
% Im
%   C:[1/2, 1/2]
%   T:[cos(d), 1]
% Ip will remain if 180-pulse is not calibrated.