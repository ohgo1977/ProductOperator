clear
close all

Ix = PO(3,{'Ix'}).M;
Iy = PO(3,{'Iy'}).M;
Iz = PO(3,{'Iz'}).M;
Sx = PO(3,{'Sx'}).M;
Sy = PO(3,{'Sy'}).M;
Sz = PO(3,{'Sz'}).M;
Kx = PO(3,{'Kx'}).M;
Ky = PO(3,{'Ky'}).M;
Kz = PO(3,{'Kz'}).M;

H = 4*Ix*Sx;
r = 4*Iy*Sz;

% H = 4*Ix*Sy;
% r = 4*Iy*Sz*Ky;


% Ix = PO(4,{'Ix'}).M;
% Iy = PO(4,{'Iy'}).M;
% Sy = PO(4,{'Sy'}).M;
% Sz = PO(4,{'Sz'}).M;
% Kx = PO(4,{'Kx'}).M;
% Lx = PO(4,{'Lx'}).M;
% H = 2*Ix*Sy;
% r = 8*Iy*Sz*Kx*Lx;