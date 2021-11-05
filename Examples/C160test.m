clear
close all

syms o1 J12 t1 Zpfg 
sp = (exp(-o1*t1*1i)*exp(-pi*J12*t1*1i)*1i)/8 + (exp(-o1*t1*1i)*exp(pi*J12*t1*1i)*1i)/8 - (exp(o1*t1*1i)*exp(-pi*J12*t1*1i)*1i)/8 - (exp(o1*t1*1i)*exp(pi*J12*t1*1i)*1i)/8 - (exp(-2*Zpfg)*exp(o1*t1*1i)*exp(-pi*J12*t1*1i)*1i)/8 - (exp(-2*Zpfg)*exp(o1*t1*1i)*exp(pi*J12*t1*1i)*1i)/8 + (exp(3.5*Zpfg)*exp(-o1*t1*1i)*exp(-pi*J12*t1*1i)*1i)/8 + (exp(2*Zpfg)*exp(-o1*t1*1i)*exp(pi*J12*t1*1i)*1i)/8;

st = char(sp);
id_Zpfg = findstr(st,'Zpfg');
id_lp = findstr(st,'(');
id_rp = findstr(st,')');
id_exp = findstr(st,'exp');

PFG_mat = [];
for jj = 1:length(id_Zpfg)
    id_tmp = id_Zpfg(jj);

    id_lp_vec = id_lp(find(id_lp < id_tmp));
    id_lp_tmp = id_lp_vec(end);

    id_rp_vec = id_rp(find(id_rp > id_tmp));
    id_rp_tmp = id_rp_vec(1);

    id_exp_vec = id_exp(find(id_exp < id_tmp));
    id_exp_tmp = id_exp_vec(end);

    cv1 = id_lp_vec(find(id_lp_vec > id_exp_tmp));
    cv2 = id_lp_vec(find(id_lp_vec < id_lp_tmp));
    Lia = ismember(cv1,cv2);

    id_rp_tmp2 = id_rp_vec(1+sum(Lia));
    st_tmp = st(id_exp_tmp : id_rp_tmp2);
    fprintf(1,'%s\n',st_tmp);
    PFG_mat = cat(2,PFG_mat,str2sym(st_tmp));
end
PFG_mat = unique(PFG_mat);

coef_out = sp;
for jj = 1:length(PFG_mat)
    coef_out = subs(coef_out,PFG_mat(jj),0);
end
coef_out = simplify(coef_out,10);