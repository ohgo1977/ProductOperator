clear
close all

syms f b
H = PO(3,{'Ix' 'Iy'},{cos(f) sin(f)});
expH1 = simplify(expm(-1i*b*H.M));

Ix = PO(3,{'Ix'});
Iy = PO(3,{'Iy'});
Iz = PO(3,{'Iz'});

expH2 = simplify(expm(-1i*f*Iz.M)*expm(-1i*b*Ix.M)*expm(1i*f*Iz.M));

simplify(expH1 - expH2)