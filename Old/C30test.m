clear
close all

syms B q
rho = PO(4,{'I1z' 'I2z' 'I3z' 'S4z'},{B B B sym(1)},{'I1' 'I2' 'I3' 'S4'});
rho = simpulse(rho,{'I*' 'S4'},{'x' 'x'},{3/2*pi pi});
rho = simjc(rho,{'I1S4' 'I2S4' 'I3S4'},{q q q});

