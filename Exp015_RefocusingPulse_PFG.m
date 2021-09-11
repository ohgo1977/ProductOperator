clear
close all
% Keeler, J., Understanding NMR Spectroscopy, p. 406, 11.12.3
% Gradient G - 180+d pulse - Gradient G 

syms G gH d
rho = PO(1,{'Ix' 'Iy'},{sym(1) sym(1i)});% Ip
rho = pfg(rho, G, {gH});
rho = simpulse(rho,{'I*'},{'x'},{pi+d});
rho = pfg(rho, G, {gH});

Ip_coef = simplify(1/2*rho.coef(1) + 1/(2*1i)*rho.coef(2));
Im_coef = simplify(1/2*rho.coef(1) - 1/(2*1i)*rho.coef(2));

rho2 = dispRaiseLower(rho);
