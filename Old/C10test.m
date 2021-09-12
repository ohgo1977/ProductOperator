clear
close all

spin_label_cell = {'I1' 'I2' 'I3' 'S4'};

%sp_cell = {'I1' 'I2' 'I3' 'S4'};
%ph_cell = {'y' 'y' 'y' 'y'};
%q_cell = {pi/2 pi/2 pi/2 pi/2};

%sp_cell = {'I*' 'S4'};
%ph_cell = {'y' 'y'};
%q_cell = {pi/2 pi/2};

%sp_cell = {'*'};
%ph_cell = {'y'};
%q_cell = {pi/2};

sp_cell = {'I*' 'S*'};
ph_cell = {'y' 'y'};
q_cell = {pi/2 pi/2};

for ii = 1:max(size(sp_cell))
   sp = sp_cell{ii};
   ph = ph_cell{ii};
   q = q_cell{ii};

   if ~isempty(strfind(sp,'*'))% 'I*' or '*'
      if length(sp) == 1% '*'
        for jj = 1:max(size(spin_label_cell))
            sp_tmp = spin_label_cell{jj};
               % obj_tmp = pulse(obj_tmp,sp_tmp,ph,q);%Run pulse each
               fprintf(1,'%s %s %g\n\n',sp_tmp,ph,q);
        end
      elseif length(sp) == 2
        id_vec = find(contains(spin_label_cell,sp(1)));
        for jj = id_vec
            sp_tmp = spin_label_cell{jj};
               % obj_tmp = pulse(obj_tmp,sp_tmp,ph,q);%Run pulse each
               fprintf(1,'%s %s %g\n\n',sp_tmp,ph,q);
        end
      end
   else % sp doesn't include '*'
   % obj_tmp = pulse(obj_tmp,sp,ph,q);%Run pulse each
   fprintf(1,'%s %s %g\n\n',sp,ph,q); 
   end
end