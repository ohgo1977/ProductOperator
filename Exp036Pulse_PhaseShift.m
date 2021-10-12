clear
close all

% spin_label_cell = {'I'};
spin_label_cell = {'I' 'S'};

PO.create(spin_label_cell);
fprintf(1,'Evolution of Iz under a pulse with flip angle q and phase f\n')
fprintf(1,'Number of Spins: %d \n',length(spin_label_cell));

rho = Iz;
rho.disp = 0;

tic;
obj1 = rho.pulse_phshift('I',f,q);% UrhoUinv_mt() is called
et1 = toc;
fprintf(1,'UrhoUinv_mt(): %g s\n',et1);

H = q*(Ix*cos(f) + Iy*sin(f));
tic;
obj2 = UrhoUinv(rho,H,1);% UrhoUinv_M() is called
et2 = toc;
fprintf(1,'UrhoUinv_M(): %g s\n',et2);
