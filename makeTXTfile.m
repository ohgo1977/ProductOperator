clear

fname_list=ls('Exp*.m');

for ii = 1:size(fname_list,1)
    fname = deblank(fname_list(ii,:));
    txt_fname = sprintf('%s.txt',fname(1:end-2));
    copyfile(fname,txt_fname);
end
