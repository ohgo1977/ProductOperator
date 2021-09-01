clear
close all
% Homonuclear INEPT
% Movellan, T.K., ..., Andreas, L. B.
% J. Am. Chem. Soc. 2020, 142, 2704-2078.

% Homonuclear pulses thus the phases of simpulse() shoud be same
phid = 1:1;
ph1tab = [2 2 0 0];% Converted from (1H, 15N) phases for CP => 15N One pulse phase
ph2tab = [0*ones(1,8) 1*ones(1,8)];
ph3tab = [0*ones(1,16) 1*ones(1,16)];
ph4tab = [0 2];
ph5tab = [1*ones(1,4) 3*ones(1,4)];
phRtab = [1 3 3 1 3 1 1 3 3 1 1 3 1 3 3 1 3 1 1 3 1 3 3 1 1 3 3 1 3 1 1 3];

% Symbolic constants
syms B J t oI oS t1

coef = [];
for ii = phid
    fprintf(1,'%2d\n',ii)
    ph1 = PO.phmod(ph1tab,ii);
    ph2 = PO.phmod(ph2tab,ii);
    ph3 = PO.phmod(ph3tab,ii);
    ph4 = PO.phmod(ph4tab,ii);
    ph5 = PO.phmod(ph5tab,ii);
    phR = PO.phmod(phRtab,ii);
    
% Short CP: only I spin being close to 1Hs is polarized.
    rho = PO(2,{'Iz'});% Both I and S are 15N
    rho.dispPOtxt();
    rho = pulse(rho,'I',ph1,pi/2);

% Long CP: both I and S spins are excited.
%     rho = PO(2,{'Iz' 'Sz'});% Both I and S are 15N
%     rho.dispPOtxt();
%     rho = simpulse(rho,{'I' 'S'},{ph1 ph1},{pi/2 pi/2});

    % 1st INEPT
    rho = jc(rho,'IS',pi*J*t);
    rho = simpulse(rho,{'I' 'S'},{ph2 ph2},{pi pi});
    rho = jc(rho,'IS',pi*J*t);

    % 90 pulse - t1 - 90 pulse
    rho = simpulse(rho,{'I' 'S'},{'y' 'y'},{pi/2 pi/2});
    rho = simcs(rho,{'I' 'S'},{oI*t1 oS*t1});
    id_vec = findcoef(rho,{sin(oI*t1) sin(oS*t1)});
    rho = delPO(rho,id_vec);% Delete the term with sin(oI*t1) and sin(oS*t1)
    rho = simpulse(rho,{'I' 'S'},{'y' 'y'},{pi/2 pi/2});
    
    % 2nd INEPT
    rho = jc(rho,'IS',pi*J*t);
    rho = simpulse(rho,{'I' 'S'},{ph3 ph3},{pi pi});
    rho = jc(rho,'IS',pi*J*t);

    % Z-filter
    rho = simpulse(rho,{'I' 'S'},{ph4 ph4},{pi/2 pi/2});
    rho = simpulse(rho,{'I' 'S'},{'x' 'x'},{pi/2 pi/2});

    rho = delPO(rho,{'IxSz'});% delete 2IxSz term
    
    % 15N => 1H CP
    % ph5 is y or -y
    % 180 phase shift of ph5 changes the sign of the signal amplitude.
    if ph5 == 1
        ph5sign = 1;
    elseif ph5 == 3
        ph5sign = -1;
    end

    % Receiver
    % phR is y or -y
    % 180 phase shift of phR changes the sign of the signal amplitude.
    if phR == 1
        phRsign = 1;
    elseif phR == 3
        phRsign = -1;
    end
    
    coefI_tmp = rho.coef(1)*ph5sign*phRsign;

    I_tmp = coeffs(coefI_tmp,cos(oI*t1));
    I_tmp = I_tmp(2);
     
    S_tmp = coeffs(coefI_tmp,cos(oS*t1));
    S_tmp = S_tmp(2);
    
    coef = cat(1,coef,simplify([I_tmp S_tmp],100));
        
%     coefS_tmp = rho.coef(1)*ph5sign*phRsign;
%     coefI_tmp = rho.coef(2)*ph5sign*phRsign;
% 
%     I_tmp = coeffs(coefS_tmp + coefI_tmp,cos(oI*t1));
%     I_tmp = I_tmp(2);
%      
%     S_tmp = coeffs(coefS_tmp + coefI_tmp,cos(oS*t1));
%     S_tmp = S_tmp(2);
%     
%     coef = cat(1,coef,[I_tmp S_tmp]);
end
