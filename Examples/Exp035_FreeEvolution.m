% Exp035_FreeEvolution.m
% Comparison of the calculation speeds between UrhoUinv_mt() and UrhoUinv_M()

clear
close all

% spin_label_cell = {'I' 'S'};% Case of two spins
spin_label_cell = {'I' 'S' 'K'};% Case of three spins
PO.create(spin_label_cell);

rho = Ix;
rho.disp = 0;
fprintf(1,'Evolution of Ix under chemical shift and J-coupling\n')
fprintf(1,'Number of Spins: %d\n',length(spin_label_cell));

tic;
obj1 = rho.cs({'I'},{o1*t}).jc({'IS'},{pi*JIS*t});% UrhoUinv_mt() is called
et1 = toc;
fprintf(1,'UrhoUinv_mt(): %g s\n',et1);
dispPOtxt(obj1);

H = o1*Iz + pi*JIS*2*Iz*Sz;
tic;
obj2 = UrhoUinv(rho,H*t,1);% UrhoUinv_M() is called
et2 = toc;
fprintf(1,'UrhoUinv_M(): %g s\n',et2);
dispPOtxt(obj2);

% Rewrite obj2.coef
coef_new = simplify(rewrite(obj2.coef,'sincos'));
obj3 = set_coef(obj2,coef_new);
fprintf(1,'Rewrite coef from UrhoUinv_M()\n')
dispPOtxt(obj3);