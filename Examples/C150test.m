clear
close all

syms G Z d gH
s = - exp(-G*Z*gH*2i)/2 + (exp(-d*2i)*exp(G*Z*gH*3i))/4 + exp(d*2i)/4;
s2 = subs(s,G*Z*gH,1i*Z);

splim = limit(s2,Z,inf);
smlim = limit(s2,Z,-inf);
s_out1 = (splim + smlim)/2;
s_out2 = subs(s_out1, Inf, 0);


% syms G Z d gH
% s = - exp(-G*Z*gH*2i)/2 + (exp(-d*2i)*exp(G*Z*gH*3i))/4 + exp(d*2i)/4;
% sp = subs(s,G*Z*gH,1i*Z);
% sm = subs(s,G*Z*gH,-1i*Z);

% splim = limit(sp,Z,inf);
% smlim = limit(sm,Z,inf);
% s_out1 = (splim + smlim)/2;
% s_out2 = subs(s_out1, Inf, 0);
