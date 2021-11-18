% ------------------------------------------------------------------------
% Class       : PO
% Description : Functions for Product Operator Formalism of spin-1/2
% Requirement : MATLAB Symbolic Math Toolbox
% Developer   : Dr. Kosuke Ohgo
% ULR         : https://github.com/ohgo1977/ProductOperator
% Version     : 1.0.0
%
% Please read the manual (PO_Manual.pdf) for details.
%
% ------------------------------------------------------------------------
%
% MIT License
%
% Copyright (c) 2021 Kosuke Ohgo
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%%
classdef (InferiorClasses = {?sym}) PO < matlab.mixin.CustomDisplay
% Overloading some functions of sym class, custum display of properties

    % https://www.mathworks.com/help/matlab/matlab_oop/property-attributes.html
    %%    
    properties (SetAccess = protected) % Read-Only from the Command Window
        axis        % Showing the status of axis direction for each spin.
                    % 0:E 1:x 2:y 3:z 4:p 5:m 6:a 7:b
                    % The column size corresponds to the number of spin types in the system.

        coef        % Coefficients of product operators (Symbolic).
                    % coef should not include the 2^(N-1) coefficient.

        spin_label  % Labels for spin1, 2, 3... stored in a cell. 

        basis       % String value to distinguish the basis-status in the calculations.
                    % 'xyz', 'pmz' or 'pol'

        logs        % Log information of obj. It can be turned off by changing logs_bin to 0.

    end

    %%
    properties (Access = protected) 
        % No access from the Command Window. Each PO object has its own value.
        bracket % Binary value to indicate cases with (a+b) or (a-b) type coefficient (1: yes, 0: no)

        SimplifySteps = 10
                % Number of steps used for simplify() in CombPO().
                % This prorperty can be changed from set_SimplifySteps().

        PFGq = [];
                % Angles used for PFG. PFGq is used in dephase().                
    end

    %% 
    properties (Dependent) % Parameters depending on other parameters
        Ncoef       % The 2^(N-1) coefficient for N spin-1/2 system (Symbolic) 
        txt         % Text output of Product Operators (String)
        M           % Matrix Form
        coherence   % Populations (diagonal) and coherences (off-diagonal) of a density operator
    end
    
    properties (Constant = true, Hidden = true) 
        % Constant throughout the methods. All PO objects have same value. Not displayed
        spin_label_cell_default = {'I' 'S' 'K' 'L' 'M'}; % Default labels for spin_label
        asterisk_bin = 0;
        % Control to the asterisk '*' between spin operators in the txt property.
        % 0 :w/o, 1: w/ '*'.
    end

    %%
    properties (Constant = true) 
        % Constant throughout the methods. All PO objects have same value.
        sqn = sym(1/2);% Spin Quantum Number
    end
    
    %%
    properties % No access limit
        disp = 1    % Control the display of the applied method on the monitor.
                    % 1: On, 2: Off
    end
    
    % Custom Display of properties
    % https://www.mathworks.com/matlabcentral/answers/6940-save-disp-output-as-a-string-cell-variable
    % https://www.mathworks.com/help/matlab/matlab_oop/custom-display-interface.html 
    methods (Access = protected)
        function footer = getFooter(obj)
            if ~isscalar(obj)
                footer = getFooter@matlab.mixin.CustomDisplay(obj);
            else
                footer = '';
                footer = sprintf('%sM:',footer);
                s_tmp = evalc('disp(obj.M)');
                footer = sprintf('%s\n%s',footer,s_tmp);
                footer = sprintf('%scoherence:',footer);
                s_tmp = evalc('disp(obj.coherence)');
                footer = sprintf('%s\n%s',footer,s_tmp);
            end
         end

        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = struct('txt', obj.txt);
                propList.spin_label = obj.spin_label;
                % use StructName.FieldName = FieldValue 
                % instead of struct(FieldName, FieldValue,...) to be a scalar struct.
                propList.basis = obj.basis;
                propList.disp = obj.disp;
                propList.logs = obj.logs;
                propList.axis = obj.axis;
                propList.coef = obj.coef;
                propList.Ncoef = obj.Ncoef;
                propList.sqn = obj.sqn;
                propgrp = matlab.mixin.util.PropertyGroup(propList); 
           end
        end
    end

    %%
    methods
        %% Ncoef_out = get.Ncoef(obj)
        function Ncoef_out = get.Ncoef(obj)
            if strcmp(obj.basis, 'xyz')% Cartesian operator basis
                Ncoef_out = sym(2.^(sum((obj.axis~=0),2)-1));
                % 2^(Ns-1) 
                % Examples 
                % IxSx   => Ns = 2 => Ncoef = 2
                % IxSxKx => Ns = 3 => Ncoef = 4
            elseif strcmp(obj.basis, 'pmz') || strcmp(obj.basis, 'pol')% Raising/Lowering  or Polarization operator bases
                Ncoef_out = sym(ones(size(obj.axis,1),1));% Ncoef = 1
            end
        end 
        % get.Ncoef
        
        %% txt_out = get.txt(obj)
        function txt_out = get.txt(obj)
            % txt_out = get.txt(obj)
            % Calculation of dependent Property txt.
            if isempty(find(obj.coef ~= sym(0),1))% If all coef values are zero.
                txt_out = '0';
            else

                syms A_dummy
                txt_out = '';
                char_pcoef_in = char(obj.coef);
                char_pcoef_in_cell = cell(length(obj.coef),1);
                char_mcoef_in = char(-1*obj.coef);
                char_mcoef_in_cell = cell(length(obj.coef),1);
                char_dummy_in = char(A_dummy*obj.coef);
                char_dummy_in_cell = cell(length(obj.coef),1);
                char_Ncoef_in = char(obj.Ncoef);
                char_Ncoef_in_cell = cell(length(obj.coef),1);

                if length(obj.coef) == 1
                    char_pcoef_in = ['; ' char_pcoef_in ';'];% Add ';' at the begining and end.
                    char_mcoef_in = ['; ' char_mcoef_in ';'];
                    char_dummy_in = ['; ' char_dummy_in ';'];
                    char_Ncoef_in = ['; ' char_Ncoef_in ';'];
                else
                    char_pcoef_in = ['; ' char_pcoef_in(2:end-1) ';'];% Remove '[' and ']' and Add ';' at the begining and end.
                    char_mcoef_in = ['; ' char_mcoef_in(2:end-1) ';'];
                    char_dummy_in = ['; ' char_dummy_in(2:end-1) ';'];
                    char_Ncoef_in = ['; ' char_Ncoef_in(2:end-1) ';'];
                end

                char_pcoef_in_id_tmp = strfind(char_pcoef_in,';');
                char_mcoef_in_id_tmp = strfind(char_mcoef_in,';');
                char_dummy_in_id_tmp = strfind(char_dummy_in,';');
                char_Ncoef_in_id_tmp = strfind(char_Ncoef_in,';');

                for ii = 1:length(obj.coef)
                    char_pcoef_in_cell{ii,1} = char_pcoef_in(char_pcoef_in_id_tmp(ii) + 2:char_pcoef_in_id_tmp(ii + 1) - 1);
                    char_mcoef_in_cell{ii,1} = char_mcoef_in(char_mcoef_in_id_tmp(ii) + 2:char_mcoef_in_id_tmp(ii + 1) - 1);
                    char_dummy_in_cell{ii,1} = char_dummy_in(char_dummy_in_id_tmp(ii) + 2:char_dummy_in_id_tmp(ii + 1) - 1);
                    char_Ncoef_in_cell{ii,1} = char_Ncoef_in(char_Ncoef_in_id_tmp(ii) + 2:char_Ncoef_in_id_tmp(ii + 1) - 1);
                end

                for ii = 1:length(obj.coef)
                    axis_tmp = obj.axis(ii,:);
                    pcoef_tmp = char_pcoef_in_cell{ii};
                    mcoef_tmp = char_mcoef_in_cell{ii};
                    dummy_tmp = char_dummy_in_cell{ii};
                    Ncoef_tmp = char_Ncoef_in_cell{ii};
                    bracket_tmp = obj.bracket(ii);

                    pt = axis2pt(obj,axis_tmp);
                    
                    % Remove '1' from single-type P.O. (N = 1 for 2^(N-1))
                    if ~strcmp(Ncoef_tmp,'1')
                        if obj.asterisk_bin == 0
                            ptc = strcat(Ncoef_tmp,pt);
                        elseif obj.asterisk_bin == 1
                            ptc = strcat(Ncoef_tmp,'*',pt);% with *
                        end
                    else
                        ptc = pt;
                    end

                    dummy_id = strfind(dummy_tmp,'A_dummy'); % Position of 'A_dummy'
                    minus_id = strfind(dummy_tmp,'-');       % Position of '-''

                    tf_v = 0; % 0: Positive Sign
                    if ~isempty(minus_id) && minus_id(1) < dummy_id % if '-' is on the lef of 'A_dummy'
                        tf_v = 1;% 1: Negative Sign
                    end

                    if tf_v
                        if ~strcmp(pcoef_tmp,'-1')% if coef = sym(-1), it is not required to display 1 as a coeffcieint.
                            if bracket_tmp == 1
                                ptc = strcat(ptc,'*','(',mcoef_tmp,')');% Add bracket for 'a - b'-type coefficient
                            else
                                ptc = strcat(ptc,'*',mcoef_tmp);
                            end
                        end
                        txt_out = [txt_out,' ','-',' ',ptc];% Add Negative sign in the text

                    else % Symbols with positive sign such as q in addtion to symbolic positive numbers.
                        if ~strcmp(pcoef_tmp,'1')% if coef = sym(1), it is not required to display 1 as a coeffcieint.
                            if bracket_tmp == 1
                                ptc = strcat(ptc,'*','(',pcoef_tmp,')');% Add bracket for 'a + b'-type coefficient
                            else
                                ptc = strcat(ptc,'*',pcoef_tmp);
                            end            
                        end

                        if ~strcmp(pcoef_tmp,'0')% if coef ~= sym(0)
                            if ii == 1% In the case of 1st term, no need to add '+'.
                                txt_out = [txt_out,ptc];% No Positive sign for 1st term
                            else
                                txt_out = [txt_out,' ','+',' ',ptc];
                            end
                        end
                    end
                end
            end
        end 
        % get.txt

        %% M_out = get.M(obj)
        function M_out = get.M(obj)
            % Create a Matrix representation.
            for ii = 1:size(obj.axis,1)
                axis_tmp = obj.axis(ii,:);

                for jj = 1:length(axis_tmp)
                    Ma = obj.axis2M(axis_tmp(jj),obj.sqn); 
                    if jj == 1
                        M_tmp = Ma;
                    else
                        M_tmp = kron(M_tmp, Ma);% Kronecker tensor product
                    end
                end
                % M_tmp is a double class
                % Mo below is sym class because coef and Ncoef are symbolic.
                Mo = M_tmp*obj.coef(ii)*obj.Ncoef(ii)*2^length(find(axis_tmp == 0));
                    % Example
                    % 2IxSz in 4-spin system (ISKL) (axis_tmp = [1 3 0 0])
                    % 2IxSz = 2IxSz*4*1/2E*1/2E = Ncoef*2^2*Ix*Sz*1/2E*1/2E
                    %
                    % IpSz in 4-spin system (ISKL) (axis_tmp = [4 3 0 0])
                    % IpSz = IpSz*4*1/2E*1/2E = Ncoef*2^2*Ip*Sz
 
                if ii == 1
                   M_out = Mo; 
                else
                   M_out = M_out + Mo;  
                end
            end
        end 
        % get.M
        
        %% coherence_out = get.coherence(obj)
        function coherence_out = get.coherence(obj)   
            % Populations and coherences of a density operator
            spin_no = size(obj.axis,2);
            coherence_out = PO.rho_box(spin_no);
            coherence_out(obj.M == 0) = 0;
        end
        % get.coherence

        %% obj = PO(spin_no,sp_cell,coef_cell,spin_label_cell)
        function obj = PO(spin_no,sp_cell,coef_cell,spin_label_cell) 
            % obj = PO(spin_no,sp_cell,coef_cell,spin_label_cell)
            % This is the class constructor, obj should not be involved in the argument.
            % https://www.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html
            %
            % spin_no: Number of spin types in the system.
            % For example, spin_no = 2 for I-S or I1-I2 system.
            %
            % sp_cell: Input of product operators as a cell array.
            % Example
            % sp_cell = {'Ix' 'IySz'}
            % Do not put 2^(N-1) coefficients (such as '2' for 2IySz)
            % Each term should have characters of spin and phase.
            %
            % coef_cell: Assignment of coefficients for sp_cell.
            % If only spin_no and sp_cell are input, sym(1) is automatically assigned for each term.
            % Example
            % sym q
            % {cos(q) sin(q)}
            % Do not include 2^(N-1) coefficients in coef_cell.
            %
            % spin_label_cell: Labels for spin types
            % If spin_label_cell is not input, {'I' 'S' 'K' 'L' 'M'} is used.
            % The size of spin_label_cell should be equal or larger than spin_no.
            % Example
            % spin_label_cell = {'I1' 'I2' 'I3' 'S4'}
            % In this case, sp_cell should be
            % sp_cell = {'I1x' 'I1yS4z'}
            % Note: there should be no overlap of strings among the members of spin_label_cell.
            % For example, Both 'I1' and 'I12' have the string 'I1'. This causes an error.
            %

            if nargin > 0 % This condition allows to create an empty PO object by obj = PO();
                axis_out = [];
                coef_out = [];
                bracket_out = [];

                if nargin <= 3
                    spin_label_cell = PO.spin_label_cell_default;
                end

                if length(spin_label_cell) < spin_no % Abort spin_label_cell is not large enough.
                    error('the size of spin_label_cell must be same as or bigger than spin_no');
                end

                spin_label_cell = spin_label_cell(1:spin_no);% Adjust the size of spin_label_cell to spin_no. 

                for ii = 1:max(size(sp_cell))
                    sp = sp_cell{ii};

                    axis_tmp = zeros(1,spin_no);
                    % Initial state of axis_temp is 1/2E for all spins.
                    % if sp is not listed in spin_label,
                    % a value in axis_tmp is kept as 0 then sp is considered as 1/2E.

                    for jj = 1:length(spin_label_cell)
                        spin_label_tmp = spin_label_cell{jj};
                        if contains(sp,spin_label_tmp)
                            id_tmp = jj;
                            phase_s = sp(strfind(sp,spin_label_tmp) + length(spin_label_tmp));

                            switch phase_s
                                case 'x', phase_id = 1;
                                case 'y', phase_id = 2;
                                case 'z', phase_id = 3;
                                case 'p', phase_id = 4;
                                case 'm', phase_id = 5;
                                case 'a', phase_id = 6;
                                case 'b', phase_id = 7;
                                otherwise, phase_id = 0;%Any unknown phase becomes 1/2E
                            end
                            axis_tmp(id_tmp) = phase_id;
                        end
                    end

                   axis_out = cat(1,axis_out,axis_tmp);

                   if nargin > 2
                        symcoef = sym(coef_cell{ii});
                   else
                       symcoef = sym(1);
                   end
                   coef_out = cat(1,coef_out,symcoef);

                   if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                       bracket_out = cat(1,bracket_out,1);
                   else
                       bracket_out = cat(1,bracket_out,0);
                   end
                end

                if isempty(find(axis_out == 4,1)) && isempty(find(axis_out == 5,1)) && isempty(find(axis_out == 6,1)) && isempty(find(axis_out == 7,1))
                    % If there are no p,m,a,b
                    basis_out = 'xyz';
                elseif isempty(find(axis_out == 1,1)) && isempty(find(axis_out == 2,1)) && isempty(find(axis_out == 6,1)) && isempty(find(axis_out == 7,1))
                    % If there are no x,y,a,b
                    basis_out = 'pmz';
                elseif isempty(find(axis_out == 1,1)) && isempty(find(axis_out == 2,1)) && isempty(find(axis_out == 3,1))
                    % If there are no x,y,z
                    basis_out = 'pol';
                else % Mixiture of Cartesian operator and Raising/Lowering operator bases
                    error('Error: Different operator bases (Cartesian, Raising/Lowering, Polarization) should not be mixed!!');
                end

                obj = PO(); % spin_label is empty at this point
                obj.axis = axis_out;% 1:x, 2:y, 3:z, 4:p, 5:m, 0: no type assgined
                obj.coef = coef_out;% Coefficient for the product operator
                obj.spin_label = spin_label_cell;
                obj.bracket = bracket_out;% 1: put bracket if coefficient is a sum-form.
                obj.basis = basis_out;

                obj = CombPO(obj);
                obj.logs = obj.txt;
            end
        end 
        % PO
        
        %% obj = CombPO(obj)
        function obj = CombPO(obj)
            % obj = CombPO(obj)
            % Combines coeffcieints of same type of terms in a PO-class object.
            % Also detects coefficients in which parentheses should be added and puts a flag in obj.bracket.
        
            axis_in = obj.axis;
            [~,IA,IC] = unique(axis_in,'rows');
            % Example
            % axis_in = [1 0 0; % Row 1
            %           0 2 3;  % Row 2
            %           1 0 0;  % Row 3
            %           3 3 0;  % Row 4
            %           0 2 3]; % Row 5
            % IA shows IDs of unique rows.
            % In this case, IA = [2 1 4]';
            %
            % IC is a tricky but this case IC = [2 1 2 3 1]'.
            % Then,
            % IC(1) = 2 means ROW 1 of axis_in corresponds to IA(2) = 1;
            % IC(2) = 1 means ROW 2 of axis_in corresponds to IA(1) = 2;
            % IC(3) = 2 means ROW 3 of axis_in corresponds to IA(2) = 1;
            % IC(4) = 3 means ROW 4 of axis_in corresponds to IA(3) = 4;
            % IC(5) = 1 means ROW 5 of axis_in corresponds to IA(1) = 2;
        
            coef_out = sym(zeros(length(IA),1));
            axis_out = zeros(length(IA),size(axis_in,2));
            bracket_out = zeros(length(IA),1);

            for ii = 1:length(IA) % 7 ms / ii loop
               IA_tmp = IA(ii);
               IC_tmp = find(IC == IC(IA_tmp));% IMPORTANT!!
                % IA_tmp is an ID of a unique row in axis_in.
                % Then find(IC == IC(IA_tmp)) searches IDs of rows in axis_in corresponding to the unique row.
                % Example
                % IA_tmp = IA(1) = 2. 
                % IC(IA_tmp) = IC(2) = 1. 
                % Then find(IC == IC(IA_tmp)) = find(IC == 1) = [2 5];
            
                axis_out(ii,:) = axis_in(IA_tmp,:);
                coef_out(ii) = sum(obj.coef(IC_tmp));     
            end

            coef_out = simplify(coef_out, 'Steps', obj.SimplifySteps);% Simplify with obj.SimplifySteps steps.

            A_dummy = sym('A_dummy');
            % Use "A" as a first charcter so that A_dummy comes before other alphabets.
            char_A_dummy = char(A_dummy);
            dummy_p_mat = A_dummy*coef_out;
            char_dummy_p_mat = char(dummy_p_mat);
            char_dummy_p_mat = [char_dummy_p_mat(1:end-1) ';]'];% Addition of ';]' at the end.

            sc_id = strfind(char_dummy_p_mat,';');% Postion of ';'
            A_dummy_id = strfind(char_dummy_p_mat,char_A_dummy); % Postion of 'A_dummmy'
            lp_pos = A_dummy_id + length(char_A_dummy) + 1;      % Postion of X of 'A_dummy*X'
            lp_st = char_dummy_p_mat(lp_pos);                    % Character at lp_pos
            lp_id = lp_pos(strfind(lp_st,'('));                  % id of ( of 'A_dummy*('

            bracket_out_id = [];
            for ii = 1:length(lp_id)
                lp_id_tmp = lp_id(ii);
                id_tmp = find(sc_id > lp_id_tmp);
                bracket_out_id = cat(1,bracket_out_id,id_tmp(1));
            end
            bracket_out(bracket_out_id) = 1;

            % Remove terms with 0 coefficients
            sym0 = sym(0);
            if isempty(find(coef_out ~= sym0,1)) % Only zero values in coef_out, then set as 0*1/2E. 
                id_vec = 1;
                axis_out = zeros(1,size(axis_in,2));% Reset axis_out for 1/2E
                bracket_out = 0;
            else
                id_vec = find(coef_out ~= sym0);
            end

            % To use spin_label and display in the input obj,
            % obj = PO() shoud not be used here.
            axis_out = axis_out(id_vec,:);
            coef_out = coef_out(id_vec,:);
            bracket_out = bracket_out(id_vec,:);
        
            axis_out(axis_out == 0) = 9;% Replace 0 to 9 for sorting purpose
            [axis_sort, id_sort] = sortrows(axis_out,'ascend');
            axis_sort(axis_sort == 9) = 0;% Replace 9 to 0
            obj.axis = axis_sort;
            obj.coef = coef_out(id_sort,:);
            obj.bracket = bracket_out(id_sort,:);
        end 
        % CombPO

        % obj = xyz2pmz(obj)
        function obj = xyz2pmz(obj)
            % obj = xyz2pmz(obj)
            % conversion from Cartesian operator basis to lowring/raising operator basis

            if ~strcmp(obj.basis,'xyz')
                error('The basis of the object should be xyz')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
                if isempty(find(axis_tmp == 1,1)) && isempty(find(axis_tmp == 2,1))% Only *z operators Iz, 2IzSz, 4IzSzKz,...
                    axis_out_tmp = axis_tmp;
                    coef_out_tmp = 1;
                else
                    % Conversion from Ix and Iy to Ip and Im
                    xn = length(find(axis_tmp == 1)); % Example: IxSyKxMz =>(Ip + Im)(Sp - Sm)(Kp + Km)Mz
                    yn = length(find(axis_tmp == 2)); % xn =2, yn = 1 => 2^(2+1) = 8 terms.
                    xyn = 2^(xn + yn);
                    axis_out_tmp = repmat(axis_tmp,xyn,1); % repmat([1 2 1 3],8,1)
                    axis_out_tmp(axis_out_tmp == 1) = 0;   % Remove x operators
                    axis_out_tmp(axis_out_tmp == 2) = 0;   % Remove y operators => [0 0 0 3; 0 0 0 3;... 0 0 0 3]
                    coef_out_tmp = ones(xyn,1);
            
                    dec = 0:xyn - 1;
                    bin_mat = de2bi(dec,(xn + yn),'left-msb');
                    bin_mat (bin_mat == 0) = 4;
                    bin_mat (bin_mat == 1) = 5;
                    % Creation of the pattern
                    % If xn + yn = 3, there are 8 terms with using p and m.
                    % The all combinations are
                    % 4 4 4 
                    % 4 4 5
                    % 4 5 4
                    % 4 5 5
                    % 5 4 4
                    % 5 4 5
                    % 5 5 4
                    % 5 5 5

                    int_count = 0;
                    for jj = 1:spin_no
                        axis_v = axis_tmp(jj);
                        if axis_v == 1 || axis_v == 2
                            int_count = int_count + 1;
                            bin_vec = bin_mat(:,int_count);
                            axis_out_tmp(:,jj) = bin_vec;

                            c_tmp = ones(size(bin_vec));
                            if axis_v == 1 % Ix = 1/2*Ip + 1/2*Im
                                c_tmp(bin_vec == 4) = 1/2;
                                c_tmp(bin_vec == 5) = 1/2;
                            elseif axis_v == 2 % Iy = 1/(2i)*Ip - 1/(2i)*Im
                                c_tmp(bin_vec == 4) = 1/(2*1i);
                                c_tmp(bin_vec == 5) = -1/(2*1i);
                            end
                            coef_out_tmp = coef_out_tmp.*c_tmp;
                        end
                    end
                end

                % Combine terms
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);% Ncoef => coef
                coef_out = [coef_out;coef_out_tmp];
            end % ii
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'pmz';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % xyz2pmz

        % obj = pmz2xyz(obj)
        function obj = pmz2xyz(obj)
            % obj = pmz2xyz(obj)
            % conversion from lowring/raising operator basis to Cartesian operator basis.

            if ~strcmp(obj.basis,'pmz')
                error('The basis of the object should be pmz')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
                if isempty(find(axis_tmp == 4,1)) && isempty(find(axis_tmp == 5,1))% Iz, 2IzSz, 4IzSzKz,...
                    axis_out_tmp = axis_tmp;
                    coef_out_tmp = 1;
                    xn = 0;
                    yn = 0;
                    zn = length(find(axis_tmp == 3));

                else
                    % Conversion from Ip and Im to Ix and Iy
                    xn = length(find(axis_tmp == 4));% p
                    yn = length(find(axis_tmp == 5));% m
                    zn = length(find(axis_tmp == 3));% z
                    xyn = 2^(xn + yn);
                    axis_out_tmp = repmat(axis_tmp,xyn,1);
                    axis_out_tmp(axis_out_tmp == 4) = 0; % Remove Ip
                    axis_out_tmp(axis_out_tmp == 5) = 0; % Remove Im
                    coef_out_tmp = ones(xyn,1);
            
                    dec = 0:xyn - 1;
                    bin_mat = de2bi(dec,(xn + yn),'left-msb');
                    bin_mat (bin_mat == 0) = 2;
                    bin_mat (bin_mat == 1) = 1;

                    int_count = 0;
                    for jj = 1:spin_no
                        axis_v = axis_tmp(jj);
                        if axis_v == 4 || axis_v == 5
                            int_count = int_count + 1;
                            bin_vec = bin_mat(:,int_count);
                            axis_out_tmp(:,jj) = bin_vec;

                            c_tmp = ones(size(bin_vec));
                            if axis_v == 4 % Ip = Ix + 1i*Iy
                                c_tmp(bin_vec == 1) = 1;
                                c_tmp(bin_vec == 2) = 1i;
                            elseif axis_v == 5 % Im = Ix - 1i*Iy
                                c_tmp(bin_vec == 1) = 1;
                                c_tmp(bin_vec == 2) = -1i;
                            end
                            coef_out_tmp = coef_out_tmp.*c_tmp;
                        end
                    end
                end
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii)*(1/2)^(xn + yn +zn - 1);
                % From 1*IpSpIz (Ncoef =1), IxSxIz IxSyIz etc. are created.
                % Then Ncoef for IxSxIz in xyz-basis is calculated as 4 automatically. 
                % To compensate 4 in the new Ncoef, it is necessary to apply (1/2)*(xn + yn + zn -1) to the new coef.
                % In this case (1/2)*(xn + yn + zn -1) = (1/2)*(2 + 0 + 1 -1) = 1/4.
                % This line is a main difference from the one in xyz2pmz().                    
                coef_out = [coef_out;coef_out_tmp];
            end
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'xyz';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % pmz2xyz

        % obj = xyz2pol(obj)
        function obj = xyz2pol(obj)
            % obj = xyz2pol(obj)
            % conversion from Cartesian operator basis to Polarization operator basis

            if ~strcmp(obj.basis,'xyz')
                error('The basis of the object should be xyz')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);

                % Conversion from Ix, Iy, Iz and E to Ip, Im, Ia and Ib
                xyzen = 2^spin_no;

                coef_out_tmp = ones(xyzen,1);
        
                dec = 0:xyzen - 1;
                bin_mat = de2bi(dec,spin_no,'left-msb');
                for jj = 1:spin_no
                    axis_v = axis_tmp(jj);
                    if axis_v == 1 || axis_v == 2
                        bin_mat (bin_mat(:, jj) == 0, jj) = 4; % p
                        bin_mat (bin_mat(:, jj) == 1, jj) = 5; % m
                    elseif axis_v == 3 || axis_v == 0
                        bin_mat (bin_mat(:, jj) == 0, jj) = 6; % a
                        bin_mat (bin_mat(:, jj) == 1, jj) = 7; % b
                    end
                end
                axis_out_tmp = bin_mat;

                for jj = 1:spin_no
                    axis_v = axis_tmp(jj);
                    bin_vec = bin_mat(:,jj);

                    c_tmp = ones(size(bin_vec));
                    if axis_v == 1 % Ix = 1/2*Ip + 1/2*Im
                        c_tmp(bin_vec == 4) = 1/2;
                        c_tmp(bin_vec == 5) = 1/2;

                    elseif axis_v == 2 % Iy = 1/(2i)*Ip - 1/(2i)*Im
                        c_tmp(bin_vec == 4) = 1/(2*1i);
                        c_tmp(bin_vec == 5) = -1/(2*1i);

                    elseif axis_v == 3 % Iz = 1/2*Ia - 1/2*Ib
                        c_tmp(bin_vec == 6) =  1/2;
                        c_tmp(bin_vec == 7) = -1/2;                      

                    elseif axis_v == 0 % hE = 1/2*Ia + 1/2*Ib
                        c_tmp(bin_vec == 6) =  1/2;
                        c_tmp(bin_vec == 7) =  1/2;
                    % E = Ia + Ib
                    %     c_tmp(bin_vec == 6) =  1;
                    %     c_tmp(bin_vec == 7) =  1;                                  
                    end
                    coef_out_tmp = coef_out_tmp.*c_tmp;
                end

                % Combine terms
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*2^(spin_no - 1);
                
                % If the lines above are
                % c_tmp(bin_vec == 6) =  1;
                % c_tmp(bin_vec == 7) =  1; 
                % coef_out_tmp should be calculated as
                % coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);
                %
                % (1/2)^spin_no*2^(spin_no - 1) = 2^-1
                % vs.
                % (1/2)^(spin_no - sum(axis == 0))*2^(sum(axis ~= 0) -1) =
                % 2^(sum(axis == 0) - spin_no + sum(axis ~= 0) -1) = 2^-1
                %
                % Example
                % [1 2 0 3 0]
                % 1/2*1/2*1/2*1/2*1/2*2^(spin_no - 1) = (1/2)^5*2^4 = 2^-1
                % vs.
                % 1/2*1/2*1*1/2*1*Ncoef_in = (1/2)^3*2^2 = 2^-1

                coef_out = [coef_out;coef_out_tmp];
            end % ii
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'pol';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % xyz2pol

        % obj = pol2xyz(obj)
        function obj = pol2xyz(obj)
            % obj = pol2xyz(obj)
            % conversion from Polarization operator basis to Cartesian operator basis.

            if ~strcmp(obj.basis,'pol')
                error('The basis of the object should be pol')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
                xyzen = 2^spin_no;

                coef_out_tmp = ones(xyzen,1);
        
                % Need to consider the case of Ia or Ip in multiple spins.
                % For example, Ia in three-spin system is [6 0 0] created by PO().
                % 0 in pol is converted to 0 in xyz.
                dec = 0:xyzen - 1;
                bin_mat = de2bi(dec,spin_no,'left-msb');
                for jj = 1:spin_no
                    axis_v = axis_tmp(jj);
                    if axis_v == 4 || axis_v == 5
                        bin_mat (bin_mat(:, jj) == 1, jj) = 2; % y
                        bin_mat (bin_mat(:, jj) == 0, jj) = 1; % x
                    elseif axis_v == 6 || axis_v == 7
                        bin_mat (bin_mat(:, jj) == 0, jj) = 3; % z
                        bin_mat (bin_mat(:, jj) == 1, jj) = 0; % hE
                    else
                        bin_mat(:,jj) = 0; % 0,1 => 0
                    end
                end

                axis_out_tmp = bin_mat;

                for jj = 1:spin_no
                    axis_v = axis_tmp(jj);
                    bin_vec = bin_mat(:,jj);

                    c_tmp = ones(size(bin_vec));
                    if axis_v == 4    % Ip = Ix + 1i*Iy
                        c_tmp(bin_vec == 1) = 1;
                        c_tmp(bin_vec == 2) = 1i;
                    elseif axis_v == 5 % Im = Ix - 1i*Iy
                        c_tmp(bin_vec == 1) = 1;
                        c_tmp(bin_vec == 2) = -1i;                        
                    elseif axis_v == 6 % Ia = hE + Iz
                        c_tmp(bin_vec == 0) =  1;                            
                        c_tmp(bin_vec == 3) =  1;
                    elseif axis_v == 7 % Ib = hE - Iz
                        c_tmp(bin_vec == 0) =  1;                            
                        c_tmp(bin_vec == 3) = -1;
                    end
                    coef_out_tmp = coef_out_tmp.*c_tmp;
                end
                
                axis_out = [axis_out;axis_out_tmp];
                % coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii)*(1/2)^(spin_no-1);
                coef_out_tmp = coef_out_tmp*coef_in(ii)*(1/2)^(spin_no-1);
                % IpSpKa creates IxSxKz, IxSxhE (= 1/2*IxSx), ....
                % With the automatic calculation of the new Ncoef, these terms become 
                % 4*IxSxKz instead of 1*IxSxkz (4:1) 
                % 2*IxSx   instead of 1/2*IxSx. (2:1/2 = 4:1)
                % To compensate the difference, it is necessary to apply (1/2)^(spin_no-1) to new coef.
                    
                coef_out = [coef_out;coef_out_tmp];
            end
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'xyz';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % pol2xyz

        % obj = pmz2pol(obj)
        function obj = pmz2pol(obj)
            % obj = pmz2pol(obj)
            % conversion from lowring/raising operator basis to Polarization operator basis

            if ~strcmp(obj.basis,'pmz')
                error('The basis of the object should be pol')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
                if isempty(find(axis_tmp == 3,1)) && isempty(find(axis_tmp == 0,1))% Not including Iz (3) and E (0)
                    axis_out_tmp = axis_tmp;
                    coef_out_tmp = 1;
                else
                    % Conversion from Iz and E to Ia and Ib
                    xn = length(find(axis_tmp == 0));
                    yn = length(find(axis_tmp == 3));
                    xyn = 2^(xn + yn);
                    axis_out_tmp = repmat(axis_tmp,xyn,1);
                    axis_out_tmp(axis_out_tmp == 0) = 0;
                    axis_out_tmp(axis_out_tmp == 3) = 0;
                    coef_out_tmp = ones(xyn,1);
            
                    dec = 0:xyn - 1;
                    bin_mat = de2bi(dec,(xn + yn),'left-msb');
                    bin_mat (bin_mat == 0) = 6;
                    bin_mat (bin_mat == 1) = 7;
                    % Creation of the pattern
                    % If xn + yn = 2, there are 4 terms with using a and b.
                    % The all combinations are
                    % 6 6 
                    % 6 7
                    % 7 6
                    % 7 7

                    int_count = 0;
                    for jj = 1:spin_no
                        axis_v = axis_tmp(jj);
                        if axis_v == 3 || axis_v == 0
                            int_count = int_count + 1;
                            bin_vec = bin_mat(:,int_count);
                            axis_out_tmp(:,jj) = bin_vec;

                            c_tmp = ones(size(bin_vec));
                            if axis_v == 3 % Iz = 1/2*Ia - 1/2*Ib
                                c_tmp(bin_vec == 6) = 1/2;% No Nceof in pmz, thus need to define
                                c_tmp(bin_vec == 7) = -1/2;
                            elseif axis_v == 0 % E = Ia + Ib
                                c_tmp(bin_vec == 6) = 1;% No Nceof in pmz, thus need to define
                                c_tmp(bin_vec == 7) = 1;
                            end
                            coef_out_tmp = coef_out_tmp.*c_tmp;
                        end
                    end
                end

                % Combine terms
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);% Ncoef => coef
                coef_out = [coef_out;coef_out_tmp];
            end % ii
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'pol';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % pmz2pol

        % obj = pol2pmz(obj)
        function obj = pol2pmz(obj)
            % obj = pol2pmz(obj)
            % conversion from Polarization operator basis to lowring/raising operator basis

            if ~strcmp(obj.basis,'pol')
                error('The basis of the object should be pol')
            end

            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
                if isempty(find(axis_tmp == 6,1)) && isempty(find(axis_tmp == 7,1))% Not including Ia (6) and Ib (7)
                    axis_out_tmp = axis_tmp;
                    coef_out_tmp = 1;
                else
                    % Conversion from Ia and Ib to Iz and E
                    xn = length(find(axis_tmp == 6));
                    yn = length(find(axis_tmp == 7));
                    xyn = 2^(xn + yn);
                    axis_out_tmp = repmat(axis_tmp,xyn,1);
                    axis_out_tmp(axis_out_tmp == 6) = 0;  
                    axis_out_tmp(axis_out_tmp == 7) = 0;  
                    coef_out_tmp = ones(xyn,1);
            
                    dec = 0:xyn - 1;
                    bin_mat = de2bi(dec,(xn + yn),'left-msb');
                    bin_mat (bin_mat == 0) = 0;
                    bin_mat (bin_mat == 1) = 3;
                    % Creation of the pattern
                    % If xn + yn = 2, there are 2 terms with using z and E.
                    % The all combinations are
                    % 0 0 
                    % 0 3
                    % 3 0
                    % 3 3

                    int_count = 0;
                    for jj = 1:spin_no
                        axis_v = axis_tmp(jj);
                        if axis_v == 6 || axis_v == 7
                            int_count = int_count + 1;
                            bin_vec = bin_mat(:,int_count);
                            axis_out_tmp(:,jj) = bin_vec;

                            c_tmp = ones(size(bin_vec));
                            if axis_v == 6 % Ia = 1/2*E + Iz
                                c_tmp(bin_vec == 0) = 1/2;% No Nceof with pmz basis, thus need to provide 1/2
                                c_tmp(bin_vec == 3) = 1;
                            elseif axis_v == 7 % Ib = 1/2*E + Iz
                                c_tmp(bin_vec == 0) = 1/2;% No Nceof with pmz basis, thus need to provide 1/2
                                c_tmp(bin_vec == 3) = -1;
                            end
                            coef_out_tmp = coef_out_tmp.*c_tmp;
                        end
                    end
                end

                % Combine terms
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);% Ncoef => coef
                coef_out = [coef_out;coef_out_tmp];
            end % ii
        
            bracket_out = [];
            for ii = 1:length(coef_out)
                symcoef = coef_out(ii);
                if contains(char(symcoef),'+')||contains(char(symcoef),'-')% 'a+b' or 'a-b'-type coefficients
                    bracket_out = cat(1,bracket_out,1);
                else
                    bracket_out = cat(1,bracket_out,0);
                end
            end
        
            obj2 = obj;
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.bracket = bracket_out;
            obj2.basis = 'pmz';
            obj = CombPO(obj2);
            obj.logs = char(obj2.logs,sprintf('%s',obj.txt));
        end
        % pol2pmz

        %% obj = set_basis(obj,basis_in)
        function obj = set_basis(obj,basis_in)
            % change the expression of obj using basis_in
            % set.basis causes infinite loop in the conversion methods.
            if ~strcmp(obj.basis,basis_in)
                if strcmp(obj.basis,'xyz') && strcmp(basis_in,'pmz')
                    obj = xyz2pmz(obj);
                elseif strcmp(obj.basis,'xyz') && strcmp(basis_in,'pol')
                    obj = xyz2pol(obj);
                elseif strcmp(obj.basis,'pmz') && strcmp(basis_in,'xyz')
                    obj = pmz2xyz(obj);
                elseif strcmp(obj.basis,'pmz') && strcmp(basis_in,'pol')
                    obj = pmz2pol(obj);
                elseif strcmp(obj.basis,'pol') && strcmp(basis_in,'xyz')
                    obj = pol2xyz(obj);
                elseif strcmp(obj.basis,'pol') && strcmp(basis_in,'pmz')
                    obj = pol2pmz(obj);
                end
            end
        end
        % set_basis

        %% obj = UrhoUinv_mt(obj,H,q)
        function obj = UrhoUinv_mt(obj,H,q)
            % obj = UrhoUinv_mt(obj,H,q)
            % Calculation of the evolution of rho under q*H based on the cyclic
            % commutations. The master table for the cyclic commutation is used.
            %
            % obj: PO class object with any basis. 
            %   H: PO class object with xyz-basis.
            %   q: angle in radian (symbolic or double)
            %
            % The size of H.axis should be [1,n], i.e., only single term.
            % The value in H.coef is multiplied to q in the calculation.
            % Example 1. H = -Iz
            % Rotation around -Z-axis with angle q is equivalent to rotation around Z-axis with -q
            % Exmaple 2. H = Iz*Sz 
            % H.axis = [3 3] and Ncoef = 2 is calculated from H.axis, indicating H.axis corresponds to 2IzSz.
            % Thus, H can be described as H = 2IzSz*1/2 and H.coef becomes 1/2. 
            % Keeler, Understanding NMR Spectroscopy, p. 155.
            
            % Master Table of Cyclic Commutation
            % if [A,B] = iC,[B,C] = iA, and [C,A] = iB (i.e. cyclic commutation)
            % then exp(-iqA)*B*exp(iqA) = B*cos(q) + C*sin(q)
            %             B
            %          x  y  z
            %          1  2  3
            %        ----------
            %       |C 
            %    x 1|  0  3 -2
            % A  y 2| -3  0  1
            %    z 3|  2 -1  0
            %
            % If C is 0, no actions
            % If C is a negative value, change the sign of PO.coef.

            % Check class of obj and H
            if ~isa(obj,'PO') || ~isa(H,'PO')
                error('Both obj and H should be the PO object!')
            end

            % H should be xyz basis
            if ~strcmp(H.basis,'xyz') 
                error('The basis of H should be xyz!')
            end

            % H should be a single term
            if size(H.axis,1) > 1
                error('H must be a single term!')
            end

            % H or obj should be a product of up to two operators
            if sum(H.axis ~= 0) > 2 && sum(obj.axis ~= 0) > 2 
                error('H or obj must be a product of up to two operators!')
            end

            basis_org = obj.basis;
            obj = set_basis(obj,'xyz');

            mt = [ 0  3 -2; 
                  -3  0  1; 
                   2 -1  0];

            % Conversion of q from double to Symbolic
            q = sym(q);
            q = q*H.coef;
            % Example: H.coef = -1/2, q_new = -1/2*q_old, sin(q_new) = sin(-1/2*q_old) = - sin(1/2*q_old)
            % If the line q = q*H.coef; above is not used, use the line sign_tmp = sign_tmp*sign(H.coef) below.
            % UrhoUinv_mt uses H.axis for the calculation.
            % H.axis = [3 3] automatically means 2*Iz*Sz.
            % Thus, if H = Iz*Sz, internally it is considered as 2*Iz*Sz*1/2

            % Calculation of new density operator evolved under a Hamiltonian
            H_axis = H.axis;

            type_mask_mat = (obj.axis.*H_axis) ~= 0;% Check how many spin types get matched, matched: 1, unmatched: 0
            axis_diff_mat = obj.axis ~= H_axis;% Check the difference of the direction of each spin type
            axis_mask_mat = type_mask_mat.*axis_diff_mat;
            axis_mask_vec = sum(axis_mask_mat,2);

            axis_masK_vec2 = axis_mask_vec;
            axis_masK_vec2(axis_masK_vec2 ~= 1) = 0; % Terms to be split to cos and sin: 1, otherwise: 0
            axis_new_tmp = zeros(size(obj.axis,1) + sum(axis_masK_vec2),size(obj.axis,2));
            coef_new_tmp = sym(zeros(size(obj.axis,1) + sum(axis_masK_vec2),1));

            % Put the position marker of sin-term after cos-term
            axis_mask_vec3 = axis_masK_vec2;
            id_tmp_vec = find(axis_mask_vec3 == 1);
            for ii = length(id_tmp_vec):-1:1
                id_tmp = id_tmp_vec(ii);
                if id_tmp ~= length(axis_mask_vec3)
                    axis_mask_vec3 = [axis_mask_vec3(1:id_tmp);2;axis_mask_vec3(id_tmp+1:end)];% 0: No-change, 1: cos term, 2: sin term
                elseif id_tmp == length(axis_mask_vec3)
                    axis_mask_vec3 = [axis_mask_vec3(1:id_tmp);2];
                end
            end

            % Non-sin terms
            axis_new_tmp(axis_mask_vec3 ~=2,:) = obj.axis;
            coef_new_tmp(axis_mask_vec3 ~=2) = obj.coef;
            
            % sin terms
            coef_new_tmp(axis_mask_vec3 ==2) = coef_new_tmp(axis_mask_vec3 == 1);

            % cos terms
            axis_cos = axis_new_tmp(axis_mask_vec3 == 1,:);
            coef_new_tmp(axis_mask_vec3 == 1) = coef_new_tmp(axis_mask_vec3 == 1)*cos(q);

            H_axis_mat = repmat(H_axis,size(axis_cos,1),1);

            mt_ver = 2;
            if mt_ver == 1
                mt_large = zeros(4,4);
                mt_large(1:3,1:3) = mt;
                mt_large(4,1:3) = 1:3;
                mt_large(1:3,4) = [1:3]';
                axis_cos4 = axis_cos;
                H_axis_mat4 = H_axis_mat;
                axis_cos4(axis_cos4 == 0) = 4;% 0 => 4
                H_axis_mat4(H_axis_mat4 == 0) = 4;% 0 => 4

                axis_sin = zeros(size(axis_cos));
                for ii = 1:size(axis_sin,1)
                    axis_sin_tmp = mt_large(H_axis_mat4(ii,:),axis_cos4(ii,:));
                    axis_sin(ii,:) = diag(axis_sin_tmp)';
                end

            elseif mt_ver == 2
                mt_large = zeros(4,4);
                mt_large(1:3,1:3) = mt;
                mt_large(4,1:3) = 1:3;
                mt_large(1:3,4) = [1:3]';
                axis_cos4 = axis_cos;
                H_axis_mat4 = H_axis_mat;
                axis_cos4(axis_cos4 == 0) = 4;% 0 => 4
                H_axis_mat4(H_axis_mat4 == 0) = 4;% 0 => 4

                axis_cos4 = reshape(axis_cos4',1,numel(axis_cos4));
                H_axis_mat4 = reshape(H_axis_mat4',1,numel(H_axis_mat4));
                axis_sin_tmp = mt_large(H_axis_mat4,axis_cos4);
                axis_sin_tmp = diag(axis_sin_tmp);% Column vector
                axis_sin = reshape(axis_sin_tmp',fliplr(size(axis_cos)));
                axis_sin = axis_sin';

            elseif mt_ver == 3
                axis_sin = zeros(size(axis_cos));
                for ii = 1:size(axis_cos,1)
                    for jj = 1:size(axis_cos,2)
                        if axis_cos(ii,jj) ~= 0 && H_axis_mat(ii,jj) ~= 0
                            axis_sin(ii,jj) = mt(H_axis_mat(ii,jj),axis_cos(ii,jj));
                        elseif axis_cos(ii,jj) == 0 && H_axis_mat(ii,jj) ~= 0
                            axis_sin(ii,jj) = H_axis_mat(ii,jj);
                        elseif axis_cos(ii,jj) ~= 0 && H_axis_mat(ii,jj) == 0
                            axis_sin(ii,jj) = axis_cos(ii,jj);
                        end
                    end
                end
            end

            axis_new_tmp(axis_mask_vec3 == 2,:) = axis_sin;            
            axis_new_tmp = abs(axis_new_tmp);

            if ~isempty(axis_sin)
                axis_sin_sign = axis_sin;
                axis_sin_sign(axis_sin_sign == 0) = 1;% Set the sign of 0 as 1
                axis_sin_sign = prod(sign(axis_sin_sign),2);
                coef_new_tmp(axis_mask_vec3 == 2) = coef_new_tmp(axis_mask_vec3 == 2).*axis_sin_sign*sin(q);
            end

            coef_new = coef_new_tmp;
            axis_new = axis_new_tmp;

            obj.axis = axis_new;
            obj.coef = coef_new;
            obj.bracket = zeros(size(obj.coef));
            obj = CombPO(obj);
            obj = set_basis(obj,basis_org);
        end 
        % UrhoUinv_mt

        %% obj = urhoUinv_M(obj, H, q)
        function obj = UrhoUinv_M(obj,H,q)
            % obj = UrhoUinv_M(obj,H,q)
            % The evolution of rho under q*H by the matrix calculation.
            % The speed is slower than UrhoUinv_mt()

            % Check class of obj and H
            if ~isa(obj,'PO') || ~isa(H,'PO')
                error('Both obj and H should be the PO object!!')
            end

            q = sym(q);
            obj_M1 = obj.M;
            H_M = H.M;
            obj_M2 = expm(-1i*H_M*q)*obj_M1*expm(1i*H_M*q);
            if strcmp(obj.basis, 'xyz')
                obj = PO.M2xyz(obj_M2, obj.spin_label);
            elseif strcmp(obj.basis, 'pmz')
                obj = PO.M2pmz(obj_M2, obj.spin_label);
            elseif strcmp(obj.basis, 'pol')
                obj = PO.M2pol(obj_M2, obj.spin_label);
            end 
        end
        % urhoUinv_M
        
        %% obj = UrhoUinv(obj, H, q)
        function obj = UrhoUinv(obj, H, q)
            % obj = UrhoUinv(obj, H, q)

            % Check class of obj and H
            if ~isa(obj,'PO') || ~isa(H,'PO')
                error('Both obj and H should be the PO object!!')
            end

            s0 = obj.logs;
            % switch to UrhoUinv_mt() if the conditions are met.
            if strcmp(H.basis,'xyz') && size(H.axis,1) == 1 && (sum(H.axis ~= 0) < 3 || sum(obj.axis ~= 0) < 3)
               % 1st: H should be 'xyz' basis.
               % 2nd: H should be a single term
               % 3rd: H or obj should be up to two products.
                obj = UrhoUinv_mt(obj,H,q);
                % fprintf(1,'UrhoUinv_mt\n');
            else
                obj = UrhoUinv_M(obj,H,q);
                % fprintf(1,'UrhoUinv_M\n');
            end
            obj.logs = char(s0,obj.txt);
        end
        % UrhoUinv

        %% obj = pulse1(obj,sp,ph,q)
        function [obj, id_sp] = pulse1(obj,sp,ph,q)
            % obj = pulse1(obj,sp,ph,q)
            % Calculation of the change of rho under a pulse.
            % obj: PO class object
            % sp: type of spin, character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % ph: phase of the pulse, character ('x','X','-x' etc.) or number (0,1,2,3)
            %     quadrature phase only.
            % q: flip angle in radian (symbolic or double)

            basis_org = obj.basis;
            s0 = obj.logs;
            obj = set_basis(obj,'xyz');

            axis_tmp = zeros(1,size(obj.axis,2));
            spin_label_cell = obj.spin_label;            
            if isa(sp,'char')
                for ii = 1:length(spin_label_cell)
                    if ~isempty(strfind(sp,spin_label_cell{ii}))
                        id_sp = ii;
                    end
                end
            elseif isa(sp,'double')
                id_sp = sp;
            end

            if strcmp(ph,'x')==1||strcmp(ph,'X')==1||(isa(ph,'double')&&ph == 0)
                phase_id = 1;
                coef_tmp = sym(1);
                ph_tmp = 'x';
            elseif strcmp(ph,'y')==1||strcmp(ph,'Y')==1||(isa(ph,'double')&&ph == 1)
                phase_id = 2;
                coef_tmp = sym(1);
                ph_tmp = 'y';
            elseif strcmp(ph,'-x')==1||strcmp(ph,'-X')==1||(isa(ph,'double')&&ph == 2)
                phase_id = 1;
                coef_tmp = sym(-1);
                ph_tmp = '-x';
            elseif strcmp(ph,'-y')==1||strcmp(ph,'-Y')==1||(isa(ph,'double')&&ph == 3)
                phase_id = 2;
                coef_tmp = sym(-1);
                ph_tmp = '-y';
            end
            axis_tmp(id_sp) = phase_id;

            H = PO();
            H.axis = axis_tmp;
            H.coef = coef_tmp;
            H.basis = 'xyz';
            H.bracket = 0;

            obj = UrhoUinv_mt(obj,H,q);
            obj = set_basis(obj,basis_org);

            if isa(q,'sym') == 1
                s_out = sprintf('Pulse: %s %s %s',spin_label_cell{id_sp},char(q),ph_tmp);
            else
                s_out = sprintf('Pulse: %s%4d%s',spin_label_cell{id_sp},round(q/pi*180),ph_tmp);
            end

            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end
        end 
        % pulse
        
        %% obj = pulse(obj,sp_cell,ph_cell,q_cell)
        function obj = pulse(obj,sp_cell,ph_cell,q_cell)
            % obj = pulse(obj,sp_cell,ph_cell,q_cell)
            % Calculation of the change of rho under simultaneous pulses.
            % obj: PO class object
            % sp_cell: type of spins in a cell.
            % Character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % ph_cell: phases of the pulses in a cell. 
            % Character ('x','X','-x' etc.) or number (0,1,2,3)
            % quadrature phase only.
            % q_cell: flip angles (rad) in a cell (symbolic or double)
            % Exmaple:
            % rho_new = pulse(rho,{'I' 'S'},{'x' 'y'},{pi/2 pi/2})
            
            spin_label_cell = obj.spin_label;
            [id_vec, ii_vec] = PO.sp2id(sp_cell,spin_label_cell);
            for ii = 1:length(id_vec)
                sp = id_vec(ii); % double
                ph = ph_cell{ii_vec(ii)};
                q  =  q_cell{ii_vec(ii)};
                obj = pulse1(obj,sp,ph,q);
            end
        end 
        % simpulse
        
        %% obj = cs1(obj,sp,q)
        function [obj, id_sp] = cs1(obj,sp,q)
            % obj = cs1(obj,sp,q)
            % Calculation of the chemical shift evolution of rho.
            % obj: PO class object
            % sp: type of spin, character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % q: flip angle (symbolic or double)

            basis_org = obj.basis;
            s0 = obj.logs;
            obj = set_basis(obj,'xyz');

            axis_tmp = zeros(1,size(obj.axis,2));
            spin_label_cell = obj.spin_label;

            if isa(sp,'char')
                for ii = 1:length(spin_label_cell)
                    if contains(sp,spin_label_cell{ii})
                        id_sp = ii;
                    end
                end
            elseif isa(sp,'double')
                id_sp = sp;
            end            
                        
            axis_tmp(id_sp) = 3;%Z direction

            H = PO();
            H.axis = axis_tmp;% 1:x, 2:y, 3:z, 0: no type assgined
            H.coef = sym(1);% Coefficient for the product operator
            H.basis = 'xyz';
            H.bracket = 0;% 1: put bracket if coefficient is a sum-form.

            obj = UrhoUinv_mt(obj,H,q);

            obj = set_basis(obj,basis_org);

            if isa(q,'sym') == 1
                s_out = sprintf('CS: %s %s',spin_label_cell{id_sp},char(q));
            else
                s_out = sprintf('CS: %s%4d',spin_label_cell{id_sp},round(q/pi*180));
            end

            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end

        end 
        % cs1
        
        %% obj = cs(obj,sp_cell,q_cell)
        function obj = cs(obj,sp_cell,q_cell)
            % obj = cs(obj,sp_cell,q_cell)

            spin_label_cell = obj.spin_label;
            [id_vec, ii_vec] = PO.sp2id(sp_cell,spin_label_cell);
            for ii = 1:length(id_vec)
                sp = id_vec(ii);% double
                q  =  q_cell{ii_vec(ii)};
                obj = cs1(obj,sp,q);
            end
        end 
        % cs        
        
        %% obj = jc1(obj,sp,q)
        function [obj, id_sp] = jc1(obj,sp,q)
            % obj = jc1(obj,sp,q)

            basis_org = obj.basis;
            s0 = obj.logs;
            obj =  set_basis(obj,'xyz');% Conversion to xyz

            axis_tmp = zeros(1,size(obj.axis,2));
            spin_label_cell = obj.spin_label;

            if isa(sp,'double')
                sp_tmp = '';
                for ii = 1:2
                    id_tmp = sp(ii);
                    axis_tmp(id_tmp) = 3;%Z direction
                    sp_tmp = [sp_tmp spin_label_cell{ii}];
                end
                id_sp = sp;
            elseif isa(sp,'char')
                id_sp = zeros(1,2);
                ii_int = 0;
                for ii = 1:length(spin_label_cell)
                    spin_label_tmp = spin_label_cell{ii};
                    if contains(sp,spin_label_tmp)
                        id_tmp = ii;
                        axis_tmp(id_tmp) = 3;%Z direction
                        ii_int = ii_int + 1;
                        id_sp(ii_int) = id_tmp;
                    end
                end
                sp_tmp = sp;
            end
            
            H = PO();
            H.axis = axis_tmp;
            H.coef = sym(1);
            H.basis = 'xyz';
            H.bracket = 0;

            obj = UrhoUinv_mt(obj,H,q);

            obj =  set_basis(obj,basis_org);% Conversion to basis_org

            if isa(q,'sym') == 1
                s_out = sprintf('JC: %s %s',sp_tmp,char(q));
            else
                s_out = sprintf('JC: %s%4d',sp_tmp,round(q/pi*180));
            end

            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end
        end 
        % jc1
       
        %% obj = jc(obj,sp_cell,q_cell)
        function obj = jc(obj,sp_cell,q_cell)
            % obj = jc(obj,sp_cell,q_cell)
            obj_tmp = obj;
            for ii = 1:max(size(sp_cell))
               sp = sp_cell{ii};
               q = q_cell{ii};
               obj_tmp = jc1(obj_tmp,sp,q);
            end
            obj = obj_tmp;
        end
        % jc

        %% dispPOtxt(obj)
        function dispPOtxt(obj)
            % dispPOtxt(obj)
            % Display obj.txt
            fprintf(1,'    %s\n',obj.txt);
        end 
        % dispPOtxt
       
        %% dispPO(obj)
        function dispPO(obj)
        % dispPO(obj)
        % Display terms and corresponding coefficients.
        % Example
        % if rho is for Ix*cos(q) + 2IySz*sin(q)
        % dispPO(rho) displays
        % 1  Ix   cos(q)
        % 2 2IySz sin(q)
                fprintf(1,'\n');
        
                for ii = 1:size(obj.axis,1)
                    axis_tmp = obj.axis(ii,:);
                    pt = axis2pt(obj,axis_tmp);
                    if obj.asterisk_bin == 0
                        if strcmp(char(obj.Ncoef(ii)),'1')
                            fprintf(1,'%4d%5s%-10s %s\n',ii,'',pt,char(obj.coef(ii)));
                        else
                            fprintf(1,'%4d%5s%-10s %s\n',ii,char(obj.Ncoef(ii)),pt,char(obj.coef(ii)));                   
                        end
                    elseif obj.asterisk_bin == 1
                        if strcmp(char(obj.Ncoef(ii)),'1')
                            fprintf(1,'%4d%6s%-20s %s\n',ii,'',pt,char(obj.coef(ii)));
                        else
                            fprintf(1,'%4d%5s%s%-20s %s\n',ii,char(obj.Ncoef(ii)),'*',pt,char(obj.coef(ii)));                   
                        end
                    end
                end
                fprintf(1,'\n');
        end 
        % dispPO
       
        %% obj = pulse_phshift1(obj,sp,ph,q)
        function obj = pulse_phshift1(obj,sp,ph,q)
            % obj = pulse_phshift1(obj,sp,ph,q)
            % Calculation of the change of rho under a pulse with arbitrary phase.
            % obj: PO class object
            % sp: type of spin, character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % ph: phase of the pulse, arbitrary phase is allowed (symbolic or double).
            % q: flip angle (symbolic or double)

            % Spin Dynamics 2nd Ed., p. 252, p. 391
            %  The rotating-frame Hamiltonian of a pulse with phase ph can be described as
            % H = wnut*(Ix*cos(ph) + Iy*sin(ph)) where wnut is a nutation frequency of the pulse.
            % Then the propagator of this pulse is
            % R(q) = exp(-i*q*(Ix*cos(ph) + Iy*sin(ph))) where q = wnut*t is a flip angle.
            % R(q) can be rewritten as
            % R(q) = Rz(ph)*Rx(q)*Rz(-ph).
            % Then rho_new after the pulse is
            % rho_new = R(q)*rho*R(-q)
            % = [exp(-i*ph*Iz)*[exp(-i*q*Ix)*[exp(i*ph*Iz)*rho*exp(-i*ph*Iz)]*exp(i*q*Ix)]*exp(-i*ph*Iz)]
            %  ----3-----     ----2-----    ----1-----         ----1-----    ----2-----    ----3-----
            % meaning
            % 1. Rotation -ph along Z axis
            % 2. Rotation   q along X axis
            % 3. Rotation  ph along Z axis.
            %
            % In the case of matrix calculation, 
            % expm(-1i*q*(Ix*cos(ph) + Iy*sin(ph)))*rho*expm(1i*q*(Ix*cos(ph) + Iy*sin(ph))).
            % is equivalent to 
            % expm(-1i*ph*Iz)*expm(-1i*q*Ix)*expm(-1i*-ph*Iz)*rho*expm(1i*ph*Iz)*expm(1i*q*Ix)*expm(1i*-ph*Iz)
            %
            % Note that
            % expm(-1i*q*(Ix*cos(ph) + Iy*sin(ph)))*rho*expm(1i*q*(Ix*cos(ph) + Iy*sin(ph))) 
            % cannot be described as
            % expm(-1i*q*Ix*cos(ph))*expm(-1i*q*Iy*cos(ph))*rho*expm(1i*q*Ix*cos(ph))*expm(1i*q*Iy*cos(ph)), 
            % because [Ix,Iy] ~= 0.

            disp_org = obj.disp;
            obj.disp = 0;
            s0 = obj.logs;
            obj = cs1(obj,sp,-ph);        % 1. Rotation -ph along Z axis
            obj = pulse1(obj,sp,'x',q);   % 2. Rotation   q along X axis
            [obj, id_sp] = cs1(obj,sp,ph);% 3. Rotation  ph along Z axis.
            obj.disp = disp_org;

            spin_label_cell = obj.spin_label;
            if isa(q,'sym') == 1
                if isa(ph,'sym') == 1
                    s_out = sprintf('Pulse: %s %s %s',spin_label_cell{id_sp},char(q),char(ph));
                else
                    s_out = sprintf('Pulse: %s %s %4d',spin_label_cell{id_sp},char(q),round(ph/pi*180));
                end
            else
                if isa(ph,'sym') == 1
                    s_out = sprintf('Pulse: %s%4d %s',spin_label_cell{id_sp},round(q/pi*180),char(ph));
                else
                    s_out = sprintf('Pulse: %s%4d%4d',spin_label_cell{id_sp},round(q/pi*180),round(ph/pi*180));
                end
            end

            s1 = sprintf('%s\n',s_out);
            s2 = sprintf('    %s\n',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s',s1);
                fprintf(1,'%s',s2);
            end
        end 
        % pulse_phshift1
        
        %% obj = pulse_phshift(obj,sp_cell,ph_cell,q_cell)
        function obj = pulse_phshift(obj,sp_cell,ph_cell,q_cell)
            % obj = pulse_phshift(obj,sp_cell,ph_cell,q_cell)

            spin_label_cell = obj.spin_label;
            [id_vec, ii_vec] = PO.sp2id(sp_cell,spin_label_cell);
            for ii = 1:length(id_vec)
                sp = id_vec(ii);% double
                ph = ph_cell{ii_vec(ii)};
                q  =  q_cell{ii_vec(ii)};
                obj = pulse_phshift1(obj,sp,ph,q);
            end

        end 
        % pulse_phshift

        %% obj = pfg(obj,G,gamma_cell)
        function obj = pfg(obj,G,gamma_cell)
            % obj = pfg(obj,G,gamma_cell)
            % applys pulse field gradient to all spins.
            % G is a strengh of the field and 
            % gamma_cell a cell array including gyromagnetic ratio of the spins.
            % Symbolic constant Z is used as a stamp to show terms affected by pfg().
            % This information is used in dephase().
            % This method was obitaned from POMA by Gunter (2006).

            syms Z
            obj_tmp = obj;
            spin_label_cell = obj.spin_label;
            id_vec = 1:size(obj.axis,2);

            for jj = id_vec
                sp_tmp = spin_label_cell{jj};
                q = G*Z*gamma_cell{jj};
                obj_tmp = cs1(obj_tmp,sp_tmp,q);
                obj_tmp.PFGq = cat(2,obj_tmp.PFGq,q);
            end
            obj = obj_tmp;
        end 
        % pfg


        %% obj = dephase(obj)
        function obj = dephase(obj)
            % obj = dephase(obj)
            % delete terms affected by pfg().
            % This method was obitaned from POMA by Gunter (2006).

            basis_org = obj.basis;
            s0 = obj.logs;

            % obj = set_basis(obj,'pmz');
            % Checking execution time (s)
            % PO.create({'I' 'S' 'K'});
            % rho = pfg((Ix + Sz + Ky),3.5*G,{gH gH gC});
            % tic;dephase(rho);toc;
            %               Conversion of basis in the code (the line above)
            % input basis   xyz pmz pol Non
            %         xyz   0.7 1.0 2.6 0.6
            %         pmz   1.2 0.5 2.6 0.4
            %         pol   2.7 2.1 1.2 0.9

            PFGq_in = obj.PFGq;
            syms Zpfg
            PFGq_out = [];
            for ii = 1:length(PFGq_in)
                try % including symbolic number in se, e.g., 3 in G*Z*gH*3. Then remove that number from se.
                    se = children(PFGq_in(ii));
                    v_tmp = double(se{end});% This gets an error if se is G*Z*gH.
                    PFGq_tmp = prod([se{1:end-1}]);% Remove number
                    PFGq_out = cat(2,PFGq_out,PFGq_tmp);
                catch % not including symbolic number in se, e.g., G*Z*gH.
                    PFGq_out = cat(2,PFGq_out,PFGq_in(ii));
                end
            end
            PFGq_out = unique(PFGq_out);

            coef_in = obj.coef;
            subs_in = expand(rewrite(coef_in,'exp'));

            for jj = 1:length(PFGq_out)
                PFGq_tmp = PFGq_out(jj);
                subs_in = subs(subs_in,PFGq_tmp,-1i*Zpfg);% subs exp(-2*G*Z*gH*1i) => exp(-2*Zpfg)
            end

            % https://www.mathworks.com/matlabcentral/answers/455911-find-child-subexpressions-of-symbolic-expression
            % Unfortunately, it is very difficult to obtain the internal structure of a symbolic expression in MATLAB.
            % As a result, it is not easy to manupulate particular terms including Zpfg, i.e., exp(n*Zpfg) in, for example,
            % sp =  (exp(3.5*Zpfg)*exp(-o1*t1*1i)*exp(-pi*J12*t1*1i)*1i)/8 + (exp(2*Zpfg)*exp(-o1*t1*1i)*exp(pi*J12*t1*1i)*1i)/8;
            % In the code below, the terms exp(n*Zpfg) are extracted from the string corresponding to sb.

            % In the future version, use of regexp should be considered.

            st = char(subs_in);

            id_Zpfg = strfind(st,'Zpfg');
            id_lp = strfind(st,'(');
            id_rp = strfind(st,')');
            id_exp = strfind(st,'exp');

            PFG_mat = [];
            while length(id_Zpfg) > 0
                id_tmp = id_Zpfg(1); % Position of Zpfg
                                                                
                id_lp_vec = id_lp(id_lp < id_tmp);    % Positions of ( on the left side of Zpfg
                id_lp_tmp = id_lp_vec(end);           % Position of the closest ( on the left side of Zpfg 
                            
                id_exp_vec = id_exp(id_exp < id_tmp); % Positions of exp on the left side of Zpfg
                id_exp_tmp = id_exp_vec(end);         % Position of the closest exp on the left side of Zpfg
                                                         % exp( (n*Zpfg)/m)
                cv1 = id_lp_vec(id_lp_vec > id_exp_tmp); %  =>( ( 
                cv2 = id_lp_vec(id_lp_vec <= id_lp_tmp); %    ( (<=
                Lia = ismember(cv1,cv2);                 %    1 1, meaning there are two lps.
                                                         % exp( (n*Zpfg)/m) 
                id_rp_vec = id_rp(id_rp > id_tmp);       %             )  )
                id_rp_tmp = id_rp_vec(sum(Lia));         %                )<= additional rp to ajudst tolal two rps.

                st_tmp = st(id_exp_tmp : id_rp_tmp);
                PFG_mat = cat(2,PFG_mat,str2sym(st_tmp)); % Conversion from char to sym.
                
                % Delete exp((n*Zpfg)/m) (st_tmp) from st to reduce the number of the loop. This step is critical.
                id_st_tmp = strfind(st,st_tmp); % id_st_tmp = [m n...];
                id_tmp = id_st_tmp;
                for kk = 1:length(st_tmp)-1
                    id_tmp = cat(2,id_tmp,id_st_tmp + kk);
                end
                id_tmp = sort(id_tmp);% id_st_tmp = [m m+1 m+2 ... m+length(st_tmp)-1 n n+1 n+2 ... n+length(st_tmp)-1 ...];
                st(id_tmp) = '';
                id_Zpfg = strfind(st,'Zpfg');
                id_lp = strfind(st,'(');
                id_rp = strfind(st,')');
                id_exp = strfind(st,'exp');
            end

            PFG_mat = unique(PFG_mat);

            coef_out = subs_in;
            for jj = 1:length(PFG_mat)
                coef_out = subs(coef_out,PFG_mat(jj),0);
            end

            obj.coef = rewrite(coef_out,'sincos');

            obj = CombPO(obj);
            % obj = set_basis(obj,basis_org);% 0.1 s

            s_out = sprintf('Dephasing of terms including %s',char(PFGq_out));
            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end
        end
        % dephase


        %% obj = receiver(obj,phR)
        function obj = receiver(obj,phR)
            % obj = receiver(obj,phR)
            % Rotation around Z-axis by -phR,
            % where phR is a receiver phase.
            % This method was obitaned from POMA by Gunter (2006).

            q = PO.ph2q(phR);
            disp_org = obj.disp;
            obj.disp = 0;
            s0 = obj.logs;
            sp_cell = obj.spin_label;

            q_cell = PO.v2cell(-q,sp_cell);% Rotation of -q around Z-axis
            obj = cs(obj,sp_cell,q_cell);
            obj.disp = disp_org;

            s_out = sprintf('Receiver with %s',PO.ph_num2str(phR));
            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end
        end
        % receiver

        %% obj = observable(obj,sp_cell)
        function obj = observable(obj,sp_cell)

            spin_label_cell = obj.spin_label;
            if nargin < 2
                sp_cell = spin_label_cell;
            end

            basis_org = obj.basis;
            s0 = obj.logs;
            if ~strcmp(obj.basis,'xyz')
                obj = set_basis(obj,'xyz');
            end

            % Selection of observable terms
            axis_in = obj.axis;
            id_in = find(sum(axis_in == 1|axis_in == 2,2) == 1);% terms with only one x or one y
            obj = selPO(obj,id_in);

            id_vec = PO.sp2id(sp_cell,spin_label_cell);% target spin types

            id_in = [];
            for ii = 1:size(obj.axis,1)
                axis_tmp = obj.axis(ii,:);
                if sum(axis_tmp(id_vec) == 1|axis_tmp(id_vec) == 2) == 1% No need to take sum?
                % if axis_tmp(id_vec) == 1 || axis_tmp(id_vec) == 2
                    id_in = cat(2,id_in,ii);
                end
            end
            obj = selPO(obj,id_in);
            obj = set_basis(obj,basis_org);

            if isa(cell2mat(sp_cell),'double')
                s_out = sprintf('Selecting Observable of %s',cell2mat(spin_label_cell(cell2mat(sp_cell))));
            else
                s_out = sprintf('Selecting Observable of %s',cell2mat(sp_cell));
            end
            s1 = sprintf('%s',s_out);
            s2 = sprintf('    %s',obj.txt);
            obj.logs = char(s0,s1,s2); 
            if obj.disp == 1
                fprintf(1,'%s\n',s1);
                fprintf(1,'%s\n',s2);
            end

        end
        % observable

        %% [a0_V,rho_V] = SigAmp(obj,sp_cell,phR)
        function [a0_V,rho_V] = SigAmp(obj,sp_cell,phR)
            % a0_V = SigAmp2(obj,sp_cell,phR)
            % Calculation of initial signal amplitudes (t=0) in the equation
            % s(t) = 2*i*(rho[-b](t) + rho[-a](t) + rho[b-](t) + rho[a-](t))*exp(-i*phrec)
            % Spin Dynamics (2nd Ed.), p.379.
            %
            % Related topics: Spin Dynamics (2nd Ed.), p.262, p. 287, p. 371, p.379, pp.608-610.
            %
            % Example
            % a0_V = SigAmp(rho,{'S'},'y')
            % a0_V = SigAmp(rho,{'I' 'S'},0)
            % a0_V = SigAmp(rho,{'I*' 'S'},0)
            % a0_V = SigAmp(rho,{1 2},0)

            spin_label_cell = obj.spin_label;
            sp_cell_tmp = cell(1,0);
            ii_int = 0;
            for ii = 1:max(size(sp_cell))
                sp = sp_cell{ii};
                if isa(sp,'double')
                    sp = spin_label_cell{sp};
                end

                if contains(sp,'*')% Wildcard 'I*' or '*'
                    if length(sp) == 1 % '*'
                        id_vec = 1:size(obj.axis,2);                        
                    elseif length(sp) == 2 % 'I*' 
                        id_vec = find(contains(spin_label_cell,sp(1)));% 1st character of sp
                    end

                    for jj = id_vec
                        ii_int = ii_int + 1;
                        sp_cell_tmp{1,ii_int} = spin_label_cell{jj};% Char
                    end
                else % sp doesn't include '*'
                    ii_int = ii_int + 1;
                    sp_cell_tmp{1,ii_int} = sp;% Char
                end
            end

            for ii = 1:length(sp_cell_tmp)
                sp_tmp = sp_cell_tmp{ii};
                sp_m = [sp_tmp 'm'];% sp_m : 'Im', 'Sm', ...

                ObsPO = PO(size(obj.axis,2),{sp_m},{1},spin_label_cell);
                if ii == 1
                    obsPO_M = ObsPO.M;% Create obsPO_M
                else
                    obsPO_M = obsPO_M + ObsPO.M;
                end
            end

            obsPO_V = reshape(obsPO_M,1,numel(obsPO_M));
            id_tmp = obsPO_V ~= sym(0);
            obsPO_V = obsPO_V(id_tmp);

            objM_V = reshape(obj.M,1,numel(obj.M));
            objM_V = objM_V(id_tmp);

            a0_V = objM_V.*obsPO_V;
            % This should be Hadamard product
            % rho.M.*(Im.M + Sm.M + ...) extracts only (-1)-quantum coherence components in rho, 
            % i.e., (Im.M + Sm.M + ...) works as a mask.
            % rho.M.*(Im.M + Sm.M + ...) is equivalent to trace(rho.M*(Ip.M + Sp.M + ...))

            a0_V = 2*1i*PO.rec_coef(phR)*a0_V;

            rho = obj.coherence;

            rho_V = reshape(rho,1,numel(rho));
            rho_V = rho_V(id_tmp);

            if obj.disp == 1
                ph_s = PO.ph_num2str(phR);
                fprintf(1,'phRec: %2s\n',ph_s);
            end

        end 
        % SigAmp
       
        %% pt = axis2pt(obj,axis_tmp)
        function pt = axis2pt(obj,axis_tmp)
            % pt = axis2pt(obj,axis_tmp)
            % axis_tmp is a row vector from rho.axis.
            % obj is necessary to get spin_label from it.
            %
            % Example:
            % axis_tmp = [1 1], pt = 'IxSx'. Note that pt doesn't include '2' for '2IxSx'.
            pt = '';
            if isempty(find(axis_tmp, 1))% axis_tmp = [0 0 ... 0]
                pt = 'E';
            else
                jj_int = 0;
                for jj = 1:length(axis_tmp)
                    axis_v = axis_tmp(jj);
                    st = obj.spin_label{jj};

                    if axis_v ~=0
                        jj_int = jj_int + 1;% Counter to determine at where the 1st * should be put.
                        if axis_v == 1
                            at = 'x';  
                        elseif axis_v == 2
                            at = 'y';
                        elseif axis_v == 3
                            at = 'z';
                        elseif axis_v == 4
                            at = 'p';
                        elseif axis_v == 5
                            at = 'm';
                        elseif axis_v == 6
                            at = 'a';
                        elseif axis_v == 7
                            at = 'b';
                        end

                        if obj.asterisk_bin == 1
                            if jj_int == 1                
                                pt = strcat(pt,st,at);% No need to put * before the 1st operator.
                                                    % Example: [0 3 2] in ISK => Sz*Ky not *Sz*Ky
                            else
                                pt = strcat(pt,'*',st,at);
                            end
                        elseif obj.asterisk_bin == 0
                            pt = strcat(pt,st,at);
                        end
                    end
                end
            end
        end 
        % axis2pt     

        %% obj = delPO(obj,id_in)
        function obj = delPO(obj,id_in)
        % obj = delPO(obj,id_in)
        % rho = delPO(rho,[1 2])
        % rho = delPO(rho,{'Ix'})
        % rho = delPO(rho,{'IzSz'})
        % rho = delPO(rho,{'Ix' 'IzSz'})
        % rho = delPO(rho,{'IxS*'})
        % rho = delPO(rho,{'I*S*' 'I*S*K*'})
        % Delete particular terms in obj.
        % id_in is a vector including id numbers for the terms to be deleted.
        % These number can be obtained by dispPO(obj).

            s0 = obj.logs;
            id_in = obj.findterm(id_in);
            obj.axis(id_in,:) = [];
            obj.coef(id_in,:) = [];
            obj = CombPO(obj);
            obj.logs = char(s0,obj.txt);
        end 
        % delPO

        %% obj = selPO(obj,id_in)
        function obj = selPO(obj,id_in)
        % obj = selPO(obj,id_in)
        %
        % rho = selPO(rho,[1 2])
        % rho = selPO(rho,{'Ix'})
        % rho = selPO(rho,{'IzSz'})
        % rho = selPO(rho,{'Ix' 'IzSz'})
        % rho = selPO(rho,{'IxS*'})
        % rho = selPO(rho,{'I*S*' 'I*S*K*'})
        % Select particular terms in obj.
        % id_in is a vector including id numbers for the terms to be seleted.
        % These number can be obtained by dispPO(obj).

            s0 = obj.logs;
            id_in = obj.findterm(id_in);
            obj.axis = obj.axis(id_in,:);
            obj.coef = obj.coef(id_in,:);
            obj = CombPO(obj);
            obj.logs = char(s0,obj.txt);
        end 
        % selPO

        %% id_out = findterm(obj,id_in)
        function id_out = findterm(obj,id_in)
            % id_out = findterm(obj,id_in)
            %
            % id_out = findterm(rho,[1 2])
            % id_out = findterm(rho,{'Ix'})
            % id_out = findterm(rho,{'IzSz'})
            % id_out = findterm(rho,{'Ix' 'IzSz'})
            % id_out = findterm(rho,{'IxS*'})
            % id_out = findterm(rho,{'I*S*' 'I*S*K*'})
            
            if isa(id_in,'cell') % case of cell
                spin_label_cell = obj.spin_label;
                spin_no = size(obj.axis,2);
                sp_cell = id_in;
                id_in = []; % Initialize id_in

                for ii = 1:max(size(sp_cell))
                    sp = sp_cell{ii};

                    axis_tmp = zeros(1,spin_no); % Empty vector 
                    for jj = 1:length(spin_label_cell)

                        spin_label_tmp = spin_label_cell{jj};
                        if contains(sp,spin_label_tmp)% search sp in spin_label_cell
                            id_tmp = jj;% Column ID of axis, i.e. each spin-type
                            phase_s = sp(strfind(sp,spin_label_tmp) + length(spin_label_tmp));% Phase character
                            switch phase_s
                                case 'x', phase_id = 1;
                                case 'y', phase_id = 2;
                                case 'z', phase_id = 3;
                                case 'p', phase_id = 4;
                                case 'm', phase_id = 5;
                                case 'a', phase_id = 6;
                                case 'b', phase_id = 7;
                                case '*', phase_id = [1:7]';% Wildcard
                                otherwise, phase_id = 0;
                            end

                            if length(phase_id) == 1
                                    axis_tmp(:,id_tmp) = phase_id;
                            elseif length(phase_id) == 7% Wildcard for phase
                                rate_tmp = size(axis_tmp,1);%
                                axis_tmp = repmat(axis_tmp,7,1);% Expand the row size of axis_tmp 5 times.
                                phase_id_vec = repelem(phase_id, rate_tmp);
                                % phase_id_vec = [1 1 1... 2 2 2... 3 3 3... ... 7 7 7...]'
                                axis_tmp(:,id_tmp) = phase_id_vec;                            
                            end
                        end
                    end % for jj

                    % Detect the IDs of obj.axis corresponding to axis_tmp
                    id_tmp2 = find(ismember(obj.axis,axis_tmp,'row'));% returns column
                    id_in = cat(1,id_in,id_tmp2);% connect columns
                end % for ii
            end % isa(id_in,'cell')
            id_out = id_in;
        end
        % findterm

        %% id_vec = findcoef(obj,coef_in_cell)
        function id_vec = findcoef(obj,coef_in_cell)
            % id_vec = findcoef(obj,coef_in_cell)
            % find terms with coefficients assigned as coef_in_cell
            % coef_in_cell is a cell including coefficients (symbolic).
            % Example:
            % coef_in_cell = {sin(oI*t1) sin(oS*t1)}

            id_vec = [];
            for jj = 1:max(size(coef_in_cell))
                coef_in = coef_in_cell{jj};
                for ii = 1:length(obj.coef)
                    if has(obj.coef(ii),coef_in)
                        id_vec = cat(2,id_vec,ii);
                    end 
                end
            end
        end 
        % findcoef

        %% obj3 = commutator(obj1, obj2)
        function obj3 = commutator(obj1, obj2)
            % obj3 = commutator(obj1, obj2)
            % Commutation between obj1 and obj2.
            % If obj1 and obj2 don't commute, obj3 is 0
            % if [A,B] = iC, then B ==> B*cos(q) + C*sin(q) under A. 

            obj3 = obj1*obj2 - obj2*obj1;
        end 
        % commutator

        %% dispProp(obj, PropertyName)
        function dispProp(obj, PropertyName)
            % dispProp(obj, PropertyName)
            % Display a property PropertyName in obj.
            % PropertyName is a string.
            % This function is useful to check a protected property.
            % There may be a MATLAB-native method to check a protected property.

                prop_out = eval(['obj.',PropertyName]);
                disp(prop_out)
        end
        % dispProp

        %% obj = set_SimplifySteps(obj,new_v)
        function obj = set_SimplifySteps(obj,new_v)
            % obj = set_SimplifySteps(obj,new_v)
            % Change the property SimplifySteps to a new value.
            obj.SimplifySteps = new_v;
        end
        % set_SimplifySteps

        %% obj = set_coef(obj,new_v)
        function obj = set_coef(obj,new_v)
            % obj = set_coef(obj,new_v)
            % Change the property coef to a new value.
            obj.coef = new_v;
            obj = CombPO(obj);
        end
        % set_coef

        %% obj = plus(obj1, obj2)
        function obj = plus(obj1, obj2)
            if isa(obj2,'double')||isa(obj2,'sym') ||isa(obj2,'char') % obj1 + a
                obj_base = obj1;
                axis_tmp = zeros(1,size(obj_base.axis,2));
                if strcmp(obj1.basis,'xyz')
                    coef_tmp = 2*sym(obj2);
                else
                    coef_tmp = sym(obj2);
                end
                bracket_tmp = 0;

            elseif isa(obj1,'double')||isa(obj1,'sym') ||isa(obj1,'char') % a + obj2
                obj_base = obj2;
                axis_tmp = zeros(1,size(obj_base.axis,2));
                if strcmp(obj2.basis,'xyz')
                    coef_tmp = 2*sym(obj1);
                else
                    coef_tmp = sym(obj1);
                end
                bracket_tmp = 0;
    
            elseif isa(obj1, 'PO') && isa(obj2, 'PO')
                if size(obj1.axis,2) ~= size(obj2.axis,2)
                    error('The number of spin types for obj1 and obj2 must be same!')
                end

                % xyz + pmz => pmz
                if (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pmz')) || (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'xyz'))
                    obj1 = set_basis(obj1,'pmz');
                    obj2 = set_basis(obj2,'pmz');

                % xyz + pol => pol
                % pmz + pol => pol
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz')) || ...
                       (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    obj1 = set_basis(obj1,'pol');
                    obj2 = set_basis(obj2,'pol');
                end

                obj_base = obj1;
                axis_tmp = obj2.axis;
                coef_tmp = obj2.coef;
                bracket_tmp = obj2.bracket;              

            end
            axis_new = [obj_base.axis; axis_tmp];
            coef_new = [obj_base.coef; coef_tmp];
            bracket_new = [obj_base.bracket; bracket_tmp];

            obj_base.axis = axis_new;
            obj_base.coef = coef_new;
            obj_base.bracket = bracket_new;
            obj = CombPO(obj_base);
            obj.logs = obj.txt;
        end
        % plus

        %% obj = minus(obj1, obj2)
        function obj = minus(obj1, obj2)
            if isa(obj2,'double') || isa(obj2,'sym') || isa(obj2,'char') % obj1 - a
                obj_base = obj1;
                axis_tmp = zeros(1,size(obj_base.axis,2));
                if strcmp(obj1.basis,'xyz')
                    coef_tmp = -2*sym(obj2);
                else
                    coef_tmp = -sym(obj2);
                end
                bracket_tmp = 0;

            elseif isa(obj1,'double') || isa(obj1,'sym') || isa(obj1,'char') % a - obj2
                obj_base = obj2;
                obj_base.coef = -1*obj_base.coef;% Difference from plus()
                axis_tmp = zeros(1,size(obj_base.axis,2));
                if strcmp(obj2.basis,'xyz')
                    coef_tmp = 2*sym(obj1);
                else
                    coef_tmp = sym(obj1);
                end
                bracket_tmp = 0;
    
            elseif isa(obj1, 'PO') && isa(obj2, 'PO')
                if size(obj1.axis,2) ~= size(obj2.axis,2)
                    error('The number of spin types for obj1 and obj2 must be same!')
                end

                % xyz + pmz => pmz
                if (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pmz')) || (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'xyz'))
                    obj1 = set_basis(obj1,'pmz');
                    obj2 = set_basis(obj2,'pmz');

                % xyz - pol => pol
                % pmz - pol => pol
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz')) || ...
                       (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    obj1 = set_basis(obj1,'pol');
                    obj2 = set_basis(obj2,'pol');

                end
                obj_base = obj1;
                axis_tmp = obj2.axis;
                coef_tmp = -1*obj2.coef;% Difference from plus()
                bracket_tmp = obj2.bracket;              

            end
            axis_new = [obj_base.axis; axis_tmp];
            coef_new = [obj_base.coef; coef_tmp];
            bracket_new = [obj_base.bracket; bracket_tmp];

            obj_base.axis = axis_new;
            obj_base.coef = coef_new;
            obj_base.bracket = bracket_new;
            obj = CombPO(obj_base);
            obj.logs = obj.txt;
        end
        % minus

        %% obj = uminus(obj) 
        function obj = uminus(obj)
            obj.coef = -1*obj.coef;
            obj.logs = obj.txt;
        end
        % uminus

        %% obj = uplus(obj) 
        function obj = uplus(obj)

        end
        % uplus

        %% obj = mtimes(obj1, obj2)
        function obj = mtimes(obj1, obj2)
            if isa(obj1,'PO') && isa(obj2,'PO')
                if size(obj1.axis,2) ~= size(obj2.axis,2)
                    error('The number of spin types for obj1 and obj2 must be same!')
                end

                % xyz * xyz
                if strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'xyz')
                    basis_org = 'xyz';
                    branch_id = 1;

                % xyz * pmz or pmz * xyz 
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pmz')) || (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'xyz'))
                    basis_org = 'pmz';
                    branch_id = 2;

                % xyz * pol or pol * xyz 
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz'))
                    basis_org = 'pol';
                    branch_id = 3;

                % pmz * pmz
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pmz'))
                    basis_org = 'pmz';
                    branch_id = 4;

                % pmz * pol or pol * pmz
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    basis_org = 'pol';
                    branch_id = 5;

                % pol * pol
                elseif (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pol'))
                    basis_org = 'pol';
                    branch_id = 6;
                end
                obj1 = set_basis(obj1,'xyz');
                obj2 = set_basis(obj2,'xyz');

                % Adjustment of the row size
                obj_base = obj1;

                axis1 = obj1.axis;
                coef1 = obj1.coef;
                bracket1 = obj1.bracket;
                Ncoef1 = obj1.Ncoef;

                axis2 = obj2.axis;
                coef2 = obj2.coef;
                bracket2 = obj2.bracket;
                Ncoef2 = obj2.Ncoef;

                row1 = size(axis1,1);
                row2 = size(axis2,1);

                axis1M = repmat(axis1,row2,1);
                coef1V = repmat(coef1,row2,1);
                bracket1V = repmat(bracket1,row2,1);
                Ncoef1V = repmat(Ncoef1,row2,1);

                axis2M = repelem(axis2,row1,1);
                coef2V = repelem(coef2,row1,1);
                bracket2V = repelem(bracket2,row1,1);
                Ncoef2V = repelem(Ncoef2,row1,1);

                % Comparison of axis1M and axis2M to detect overlapping of spin types
                comp_M = axis1M.*axis2M > 0;
                comp_V = sum(comp_M,2);
                % if a value of comp_V is larger than 0, 
                % it is necessary to consider the case of IxIy, SySz,... etc.
                % axis1M         axis2M         axis1M.*axis2M             comp_V axis1M + axis2M
                % [a1 a2 a3]     [c1 c2 c3]     [1 0 0].*[0 2 0] = [0 0 0]    0   [1 0 0] + [0 2 0] = [1 2 0]
                % [b1 b2 b3]     [c1 c2 c3]      
                % [a1 a2 a3]     [d1 d2 d3]     [3 3 0].*[1 0 0] = [3 0 0]   ~0   [3 3 0] + [1 0 0] = [4 3 0]  
                % [b1 b2 b3]     [d1 d2 d3]

                % Multiplication (*) in the spin-operator system is equivalent to
                %       Addition (+) in the axis representation in this code
                % except for the multiplication/addition of the same spin type. 
                % It is because the axis value for E is 0 in the code and
                % E is the identity element of multiplication and
                % 0 is the identity element of addition.
                % Exmaple
                % Iz      : [3 0 0]
                % Sy      : [0 2 0]
                % Kx      : [0 0 1]
                % Iz*Sy*Kz: [3 2 1]

                % Iz         : [3 0 0]
                % IxSy       : [1 2 0]
                % Kx         : [0 0 1]
                % Iz*Ix*Sy*Kz: [4 2 1]
                % Need a special calculation for Iz*Ix

                axis_new = axis1M + axis2M;
                coef_new = coef1V.*coef2V;
                bracket_new = double((bracket1V + bracket2V) > 0);
                Ncoef_new = Ncoef1V.*Ncoef2V;

                at = [ 0  3  2; 
                       3  0  1; 
                       2  1  0];

                ct = [ 1/2  1i -1i;
                      -1i  1/2   1i;
                       1i -1i  1/2];
                ct = sym(1/2*ct);
                % at.*ct corresponds to
                % Ix*Ix =   1/4*E  Ix*Iy =  1i/2*Iz Ix*Iz = -1i/2*Iy
                % Iy*Ix = -1i/2*Iz Iy*Iy =   1/4*E  Iy*Iz =  1i/2*Ix
                % Iz*Ix =  1i/2*Iy Iz*Iy = -1i/2*Ix Iz*Iz =   1/4*E

                ii_vec = find(comp_V > 0);
                if ~isrow(ii_vec)
                    ii_vec = ii_vec';
                end

                for ii = ii_vec
                    comp_M_row = comp_M(ii,:);
                    jj_vec = find(comp_M_row == 1);
                    if ~isrow(jj_vec)
                        jj_vec = jj_vec';
                    end

                    for jj = jj_vec
                        id1 = axis1M(ii,jj);
                        id2 = axis2M(ii,jj);
                        axis_new(ii,jj) = at(id1,id2);
                        coef_new(ii) = coef_new(ii)*ct(id1,id2);
                    end
                end

                obj_base.axis = axis_new;
                obj_base.coef = coef_new;
                obj_base.bracket = bracket_new;

                % Comparison of old and new Ncoef
                Ncoef_in =  Ncoef_new;
                Ncoef_out = obj_base.Ncoef;
                Ncoef_rate = Ncoef_in./Ncoef_out;

                % id_coef = find(totalcoef_rate ~= sym(1));
                id_coef = find(Ncoef_rate ~= sym(1));
                if ~isrow(id_coef)
                    id_coef = id_coef';
                end

                for ii = id_coef
                    obj_base.coef(ii) = obj_base.coef(ii)*Ncoef_rate(ii);
                end

                obj = CombPO(obj_base);

                obj = set_basis(obj,basis_org);

            else
                if numel(obj1) == 1 && numel(obj2) == 1 % obj1*a or a*obj2, a should not be a vector or matrix
                    if isa(obj2,'double')||isa(obj2,'sym') ||isa(obj2,'char') % obj1*a
                        obj_base = obj1;
                        coef_tmp = sym(obj2);

                    elseif isa(obj1,'double')||isa(obj1,'sym') ||isa(obj1,'char') % a*obj2
                        obj_base = obj2;
                        coef_tmp = sym(obj1);
                    end
                    coef_new = obj_base.coef*coef_tmp;
                    obj_base.coef = coef_new;
                    obj = CombPO(obj_base);
                else % obj1*col_v or row_v*obj2
                    if (isa(obj2,'double')||isa(obj2,'sym')) && iscolumn(obj2) && size(obj2,1) == size(obj1.M,2) % obj1*col_v
                        obj = obj1.M*obj2;

                    elseif (isa(obj1,'double')||isa(obj1,'sym')) && isrow(obj1) && size(obj1,2) == size(obj2.M,1) % row_v*obj2
                        obj = obj1*obj2.M;

                    else
                        error('The sizes of obj.M and the vector should be same!');
                    end
                end
            end
            obj.logs = obj.txt;
        end
        % mtimes
        

        %% obj = mrdivide(obj1, obj2)
        function obj = mrdivide(obj1, obj2)
            if isa(obj2,'double')||isa(obj2,'sym') ||isa(obj2,'char') % obj1/a
                obj_base = obj1;
                coef_tmp = sym(obj2);
            else          
                error('PO-class object cannot the be the divisor!')
            end
            coef_new = obj_base.coef/coef_tmp;
            obj_base.coef = coef_new;
            obj = CombPO(obj_base);
            obj.logs = obj.txt;
        end
        % mrdivide
        
        %% obj = mpower(obj1, obj2)
        function obj = mpower(obj1, obj2)
            if isa(obj2,'double') && numel(obj2) == 1
                if mod(obj2,1) == 0
                    if obj2 == 0
                        obj1.coef = sym(zeros(size(obj1.coef)));
                        obj = CombPO(obj1);
                    elseif obj2 > 0
                        for ii = 1:obj2
                            if ii == 1
                                obj = obj1;
                            else
                                obj = obj*obj1;
                            end
                        end
                    else
                        error('Can''t calculate obj^-n!')
                    end
                else
                    error('n in obj^n should be rounded!')
                end
            else
                error('n in obj^n should be a scalar!')
            end
            obj.logs = obj.txt;
        end
        % mpower        
    end 
    % methods
    
    %% Static
    methods (Static)
        %% phout = phmod(phx,ii)
        function phout = phmod(phx,ii)
        % phout = phmod(phx,ii)
        % Read a phase value, phout, from a phase table, phx.
        % if ii is smaller or equal to length(phx), phout = phx(ii). 
        % If ii is larger than length(phx), then phout = phx(mod(ii,length(phx)).
        %
        % Example
        % A short phase table : phS = [0 1 2 3]
        % A  long phase table : phL = [0 1 2 3 2 3 0 1]
        % This case, there are 8 steps for complete phase cycling.
        % To repeat phS in the 8-step phase cling, phmod is used.
        % For example, phmod(phS, 5) returns a value corresponding to phS(1)

            phid = mod(ii,length(phx));
            if phid == 0
                phid = length(phx);
            end
            % phid is [1,2,3, ... length(phx)]
            phout = phx(phid);
        end 
        % phmod
       
        %% ph_s = ph_num2str(ph_n)
        function ph_s = ph_num2str(ph_n)
            % ph_s = ph_num2str(ph_n)
            % Conversion from quadrature phase number (ph_n = 0,1,2,3)
            % to phase characters (ph_s = x,y,-x,-y)

            if isa(ph_n,'char')
                ph_s = ph_n;
            else
                if ph_n == 0
                    ph_s = 'x';
                elseif ph_n == 1
                    ph_s = 'y';
                elseif ph_n == 2
                    ph_s = '-x';
                elseif ph_n == 3
                    ph_s = '-y';
                end
            end
        end 
        % ph_num2st
       
        %% coef = rec_coef(ph)
        function coef = rec_coef(ph)
            % coef = rec_coef(ph)
            % ph is a quadrature phase ('x','X',0,'y','Y',1,...)
            % then coef = exp(-1i*ph) = cos(ph) -1i*sin(ph)
            % coef(0)      =  1
            % coef(pi/2)   = -1i
            % coef(pi)     = -1
            % coef(pi*3/2) =  1i
            % Spin Dynamics (2nd Ed.), p.287.
            if strcmp(ph,'x')==1||strcmp(ph,'X')==1||(isa(ph,'double')&&ph == 0)
                coef = 1;
            elseif strcmp(ph,'y')==1||strcmp(ph,'Y')==1||(isa(ph,'double')&&ph == 1)
                coef = -1i;
            elseif strcmp(ph,'-x')==1||strcmp(ph,'-X')==1||(isa(ph,'double')&&ph == 2)
                coef = -1;
            elseif strcmp(ph,'-y')==1||strcmp(ph,'-Y')==1||(isa(ph,'double')&&ph == 3)
                coef = 1i;
            end 
        end 
        % rec_coef

        %% q = ph2q(ph)
        function q = ph2q(ph)
            % Conversion from quadrature phase to radian.
            %  x,  X, 0 => 0
            %  y,  Y, 1 => pi/2
            % -x, -X, 2 => pi
            % -y, -Y, 3 => pi*3/2

            if strcmp(ph,'x')==1||strcmp(ph,'X')==1||(isa(ph,'double')&&ph == 0)
                q = 0;
            elseif strcmp(ph,'y')==1||strcmp(ph,'Y')==1||(isa(ph,'double')&&ph == 1)
                q = pi/2;
            elseif strcmp(ph,'-x')==1||strcmp(ph,'-X')==1||(isa(ph,'double')&&ph == 2)
                q = pi;
            elseif strcmp(ph,'-y')==1||strcmp(ph,'-Y')==1||(isa(ph,'double')&&ph == 3)
                q = 3/2*pi;
            end 
        end 
        % ph2q              
       
        %% id_vec = sp2id(sp_cell,spin_label_cell)
        function [id_vec, ii_vec] = sp2id(sp_cell,spin_label_cell)
            % [id_vec, ii_vec] = sp2id(sp_cell,spin_label_cell)
            % Examples for the use of wildcard character:
            % If spin_label_cell is
            % spin_label_cell = {'I1' 'I2' 'I3' 'S4' 'S5'}
            %
            % sp_cell = {'I*' 'S*'}
            % id_vec = [1 2 3 4 5];% id of spin_label_cell
            % ii_vec = [1 1 1 2 2];% id of sp_cell
            %
            % sp_cell = {'*'}
            % id_vec = [1 2 3 4 5];% id of spin_label_cell
            % ii_vec = [1 1 1 1 1];% id of sp_cell
            %
            % sp_cell = {'I1' 'I2' 'S*'}
            % id_vec = [1 2 4 5];% id of spin_label_cell
            % ii_vec = [1 2 3 3];% id of sp_cell

            id_vec = [];
            ii_vec = [];
            for ii = 1:max(size(sp_cell))% for each element of sp_cell
                sp = sp_cell{ii};
                if isa(sp,'double')
                    id_vec_tmp = sp;

                elseif isa(sp,'char')
                    if contains(sp,'*')% Wildcard 'I*' or '*'
                        if length(sp) == 1 % '*'
                            id_vec_tmp = 1:size(spin_label_cell,2);
                        elseif length(sp) == 2 % 'I*' 
                            id_vec_tmp = find(contains(spin_label_cell,sp(1)));% 1st character of sp
                        end
                    else % sp doesn't include '*'
                        id_vec_tmp = find(contains(spin_label_cell,sp));
                    end
                end
                id_vec = cat(2,id_vec,id_vec_tmp);
                ii_vec_tmp = ii*ones(size(id_vec_tmp));
                ii_vec = cat(2,ii_vec,ii_vec_tmp);
            end
            [id_vec,id_tmp] = unique(id_vec);% Remove if duplicates exist
            ii_vec = ii_vec(id_tmp);
        end
        % sp2id

        %% M_out = axis2M(axis_v, sqn)
        function M_out = axis2M(axis_v, sqn)
            % M_out = axis2M(axis_v, sqn)
            % axis_v is a value for PO.axis, i.e., 0, 1, 2, 3, 4, 5, 6 and 7.
            % sqn is spin quantum number (Symbolic)
            % Currently, only sqn = sym(1/2) works.
            % M_out is 1/2*Pauli matrix (Symbolic)
            % (axis_v = 0 for s0 (E), 1 for s1(sx), 2 for s2(sy) and 3 for s3(sz))

            if nargin == 1
                sqn = PO.sqn;
            end
            
            if sqn == sym(1/2)
                if axis_v == 0% 1/2 E
                    M_tmp = 1/2*[1 0;0 1];
                elseif axis_v == 1% x
                    M_tmp = 1/2*[0 1;1 0];
                elseif axis_v == 2% y
                    M_tmp = 1/(2*1i)*[0 1;-1 0];
                elseif axis_v == 3% z
                    M_tmp = 1/2*[1 0;0 -1];
                elseif axis_v == 4% p
                    M_tmp = [0 1;0 0];
                elseif axis_v == 5% m
                    M_tmp = [0 0;1 0];
                elseif axis_v == 6% a
                    M_tmp = [1 0;0 0];
                elseif axis_v == 7% b
                    M_tmp = [0 0;0 1];
                end
            end
            M_out = M_tmp; % Output as double. It is converted to sym in get.M().
        end
        % axis2M

        %% rho = rho_box(n)
        function [rho, rho_num_cell] = rho_box(n)
            % Calculation of "box notation" of the density oeprator, rho.
            % a: alpha state
            % b: beta state
            % m: a => b coherence
            % p: b => a coherence
            % Spin Dynamics p. 470, p. 260, p. 160
            n_s = 2^n;
            dec = n_s-1:-1:0;
        
            bin_mat = de2bi(dec,n,'left-msb');
            % bin_mat = dec2bin(dec,n) - '0';% Alternative expression
            % Example: 3 spin system
            % a: 1, b: 0
            % |r1> = |c1> = a a a = 1 1 1
            % |r2> = |c2> = a a b = 1 1 0
            % |r3> = |c3> = a b a = 1 0 1
            % |r4> = |c4> = a b b = 1 0 0
            % |r5> = |c5> = b a a = 0 1 1
            % |r6> = |c6> = b a b = 0 1 0
            % |r7> = |c7> = b b a = 0 0 1
            % |r8> = |c8> = b b b = 0 0 0
        
            rho_cell = cell(n_s,n_s);
            % rho(c,r) = <c|rho|r>
            % Since the spin state should be read from right to left,
            % rho(c,r) correspond to |r> => |c>
            % Spin Dynamics, p. 470.
            rho_num_cell = cell(n_s,n_s);

            for ii = 1:n_s% row
                for jj = 1:n_s% column
                    r_vec = bin_mat(jj,:);%|r>
                    c_vec = bin_mat(ii,:);%|c>
                    rho_vec = c_vec - r_vec;% To know |r> => |c>
                    % rho_vec take -1, 0, 1.
                    % -1: m and 1:p
                    % 0: no change from |r> to |c>  

                    rho_tmp = char(double('b')*ones(1,n));% 'bbbb...'
                    rho_tmp(r_vec == 1) = 'a';            % 'abaa...' correspond to r_vec
                    rho_tmp(rho_vec == 1) = 'p';          % 'pbaa...'
                    rho_tmp(rho_vec == -1) = 'm';         % 'pbam...'
                    rho_cell{ii,jj} = rho_tmp;

                    rho_num_vec = zeros(size(rho_vec));
                    rho_num_vec(r_vec == 1) = 6;% a
                    rho_num_vec(rho_vec == 1) = 4;% p
                    rho_num_vec(rho_vec == -1) = 5;% m
                    rho_num_vec(rho_num_vec == 0) = 7;% b
                    rho_num_cell{ii,jj} = rho_num_vec;
                end
            end
            rho = sym(rho_cell);
        end
        % rho = rho_box(n)

        %% create(spin_label_cell)
        function create(spin_label_cell,add_cell,symcoef_switch)
            %% Spin Operators
            spin_no = length(spin_label_cell);
            assign_method = 'caller';
            for ii = 1:spin_no
                for jj = 1:7
                    if jj == 1
                        sp = [spin_label_cell{ii} 'x'];
                    elseif jj == 2
                        sp = [spin_label_cell{ii} 'y'];
                    elseif jj == 3
                        sp = [spin_label_cell{ii} 'z'];
                    elseif jj == 4
                        sp = [spin_label_cell{ii} 'p'];
                    elseif jj == 5
                        sp = [spin_label_cell{ii} 'm'];
                    elseif jj == 6
                        sp = [spin_label_cell{ii} 'a'];
                    elseif jj == 7
                        sp = [spin_label_cell{ii} 'b'];
                    end
                    obj = PO(spin_no,{sp},{1},spin_label_cell);
                    assignin(assign_method,sp,obj);
                end
            end
            obj = PO(spin_no,{'1'},{1},spin_label_cell);
            assignin(assign_method,'hE',obj);

            if nargin == 1
                symcoef_switch = 'on';
                add_cell = {};
            elseif nargin == 2
                symcoef_switch = 'on';
            end
            
            switch symcoef_switch
                case {'y','on'}
                    PO.symcoef(spin_label_cell,add_cell,'base')
                    % symcoefs stored in the workspace
                    % if PO.create is called in a function, 
                    % PO.symcoef should be also called in that function.
            end
        end
        % create

        %% symcoef(spin_label_cell)
        function symcoef(spin_label_cell,add_cell,assign_method)
            % Create pre-set Symbolic Coefficients

            spin_no = length(spin_label_cell);
            if nargin < 3
                assign_method = 'caller';
            end

            % frequency o
            for ii = 1:spin_no
                sp = spin_label_cell{ii};

                varname = ['o' sp];
                assignin(assign_method,varname,sym(varname));

                varname = ['o' sp(1)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end

                varname = ['o' sp(end)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end

                varname = ['o' num2str(ii)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end
            end

            % coupling J
            for ii = 1:spin_no
                for jj = 1:spin_no
                    if ii < jj
                        sp_ii = spin_label_cell{ii};
                        sp_jj = spin_label_cell{jj};

                        varname = ['J' sp_ii sp_jj];
                        assignin(assign_method,varname,sym(varname));

                        varname = ['J' sp_ii(1) sp_jj(1)];
                        if exist(varname,'var') == 0
                            assignin(assign_method,varname,sym(varname));
                        end

                        varname = ['J' sp_ii(end) sp_jj(end)];
                        if exist(varname,'var') == 0
                            assignin(assign_method,varname,sym(varname));
                        end

                        varname = ['J' num2str(ii) num2str(jj)];
                        if exist(varname,'var') == 0
                            assignin(assign_method,varname,sym(varname));
                        end
                    end
                end
            end
            
            % gyromagnetic ratio g
            for ii = 1:spin_no
                sp = spin_label_cell{ii};

                varname = ['g' sp];
                assignin(assign_method,varname,sym(varname));

                varname = ['g' sp(1)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end

                varname = ['g' sp(end)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end

                varname = ['g' num2str(ii)];
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end
            end

            symcoef_list = {'a' 'b' 'c' 'd' 'f' 'gH' 'gC' 'gN' 'q' 't1' 't2' 't3' 't4' 't' 'w' 'B' 'J' 'G'};
            if nargin == 1
                add_cell = {};
            end
            symcoef_list = [symcoef_list add_cell];
            for ii = 1:length(symcoef_list)
                varname = symcoef_list{ii};
                if exist(varname,'var') == 0
                    assignin(assign_method,varname,sym(varname));
                end
            end
        end
        % symcoef

        %% v_cell = v2cell(v,ref_cell)
        function v_cell = v2cell(v,ref_cell)
            % v_cell = v2cell(v,ref_cell)
            % crealt a cell array {v v v ... v} with the same size of ref_cell.
            v_cell = num2cell(v*ones(size(ref_cell)));
        end
        % v2cell

        %% line_mat = read_ascii(fname)
        function line_mat = read_ascii(fname)
            % line_mat = read_ascii(fname)
            % Read ascii file fname and output char array.

            fid=fopen(fname);
            line_mat='';
                while 1
                    tline = fgetl(fid);
                    if ~isstr(tline), break, end
                        if length(tline) == 0
                            tline=' ';
                        end
                        line_mat=strvcat(line_mat,tline);
                end
            fclose(fid);
        end
        % read_ascii

        %% [para_lines, ps_lines, rho_ini_line] = load_PS(fname)
        function [para_lines, ps_lines, rho_ini_line] = load_PS(fname)
            % [para_lines, ps_lines, rho_ini_line] = load_PS(fname)
            % Load information of pulse seuequence from fname.
            % para_lines: lines for parameters
            % ps_lines: lines for pulse sequence
            % rho_ini_line: line for rho_ini
            
            mlines = PO.read_ascii(fname);
            para_bin = 0;
            ps_bin = 0;
            rho_ini_bin = 0;
            para_lines = '';
            ps_lines = '';
            for ii = 1:size(mlines,1)
                mline = mlines(ii,:);
                if contains(mline,'% Para begin %')
                    para_bin = 1;
                    continue;
                elseif contains(mline,'% Para end %')
                    para_bin = 0;
                elseif contains(mline,'% PS begin %')
                    ps_bin = 1;
                    continue;
                elseif contains(mline,'% PS end %')
                    ps_bin = 0;
                elseif contains(mline,'rho_ini') && strfind(mline,'rho_ini') == 1
                    rho_ini_bin = 1;
                end

                if para_bin == 1
                    if rho_ini_bin == 0
                        para_lines = strvcat(para_lines,mline);
                    elseif rho_ini_bin == 1
                        rho_ini_line = mline;
                        rho_ini_bin = 0;
                    end
                elseif ps_bin == 1
                    ps_lines = strvcat(ps_lines,mline);
                end
            end
        end
        % load_PS

        %% run_PS(fname)
        function [rho_cell, rho_detect_cell, rho_total, rho_obs, a0_M, rho_M] = run_PS(fname)
            % [rho_cell, rho_total, rho_obs, a0_M, rho_M] = run_PS(fname)
            % Run a pluse sequence defined by the input file, fname.
            % The input file is an ascii file.

            %% Load lines for parameters and pulse sequence
            [para_lines, ps_lines, rho_ini_line] = PO.load_PS(fname);

            % Input parameters
            % spin_label_cell : Cell array for the spin labels. If it is not assigned, {'I' 'S' 'K' 'L' 'M'} is used.
            % rho_ini         : PO object of the initial state. It should be described by the PO objects.
            % obs_cell        : Cell array defining the observed spins. The wildcard character '*' can be used.
            %                   If it is not assigned, {'*'} is used, meaning all spins are observed.
            % ph_cell         : Cell array defining phase tables. ph_cell{n} is a row vector corresponding to phn in the pulse sequence.
            % phRtab          : Row vector defining receiver phase. Currently, only quadrature phase is accepted.
            % phid            : Row vector defining phase cycle. For full phase cycling, it is not necessary to assign this parameter.
            %                   If it is necessary to run particular sets of phases, assign this parameter. 
            %                   For example, if you like to check 2nd and 4th steps of the phase cycling, phid = [2 4];
            % coef_cell        : Cell array to create additional symbolic coefficients required for the pulse sequence. 
            %                   If it is not assigned, {} is used.
            % disp_bin        : Binary value to control the display of the pulse sequence on the Command Window (1:ON, 0:OFF).
            %                   If it is not assigned, 1 is used.

            %% Create parameters
            for ii = 1:size(para_lines,1)
                eval(para_lines(ii,:))
            end

            if exist('phid','var') == 0
                ph_length = 0;
                for ii = 1:length(ph_cell)
                    ph_length = max(ph_length,length(ph_cell{ii}));
                end
                h_length = max(ph_length,length(phRtab));
                phid = 1:ph_length;
            end

            if exist('coef_cell','var') == 0
                coef_cell = {};
            end
            
            if exist('spin_label_cell','var') == 0
                spin_label_cell = PO.spin_label_cell_default;
            end

            if exist('disp_bin','var') == 0
                disp_bin = 1;
            end

            if exist('obs_cell','var') == 0
                obs_cell = {'*'};
            end

            %% Create spin operators and sym constants
            PO.create(spin_label_cell,coef_cell,'off');
            PO.symcoef(spin_label_cell,coef_cell);% Need to call here, not in PO.create.

            %% Create rho_ini
            eval(rho_ini_line);

            %% Control of display
            rho_ini.disp = disp_bin;

            a0_M = [];
            rho_M = [];
            rho_total = 0;
            rho_cell = cell(1,length(phid));
            rho_detect_cell = cell(1,length(phid));
            int_ii = 0;
            for ii = phid
                int_ii = int_ii + 1;
                fprintf(1,'\nii: %2d\n',ii);

                % Create phases
                if exist('ph_cell','var') == 1
                    for jj = 1:length(ph_cell)
                        eval(['ph',num2str(jj),'= PO.phmod(ph_cell{jj},ii);'])
                    end
                end
                phR = PO.phmod(phRtab,ii);

                rho = rho_ini;
                if rho_ini.disp == 1
                    rho.dispPOtxt();
                end

                % Run pulse sequence
                for jj = 1:size(ps_lines,1)
                    try
                        eval(ps_lines(jj,:));
                    catch ME
                        % error_message = sprintf('PS Line %d, %s',jj,ME.message);
                        error('PS Line %d, %s',jj,ME.message);
                    end
                end

                % Store rho
                rho_cell{int_ii} = rho;

                % Receiver
                rho_detect = receiver(rho,phR);
                rho_detect_cell{int_ii} = rho_detect;
                rho_total = rho_detect + rho_total;

                if nargout > 4
                    [a0_V, rho_V] = rho.SigAmp(obs_cell,phR);
                    if nargout == 5
                        a0_M = cat(1,a0_M,a0_V);
                    end

                    if nargout == 6
                        rho_M = cat(1,rho_M,rho_V)
                    end
                end
            end
            % Observable
            rho_obs = observable(rho_total,obs_cell);
        end
        % run_PS

        %% obj = M2pol(M_in,spin_label_cell)
        function obj = M2pol(M_in,spin_label_cell)

            if ~isa(M_in,'double') && ~isa(M_in,'sym')
                error('M_in should be double or sym !')
            end

            if size(M_in,1) ~= size(M_in,2)
                error('M_in should be 2^n x 2^n !')
            end

            spin_no = log2(size(M_in,1));
            if mod(spin_no,1) ~= 0
                error('M_in should be 2^n x 2^n !')
            end

            if nargin < 2
                spin_label_cell = PO.spin_label_cell_default; % Default spin_label
            end

            if length(spin_label_cell) < spin_no % Abort spin_label_cell is not large enough.
                error('the size of spin_label_cell must be same as or bigger than spin_no');
            end

            spin_label_cell = spin_label_cell(1:spin_no);% Adjust the size of spin_label_cell to spin_no. 


            [~,rho_num_cell] = PO.rho_box(spin_no);
            M_in = sym(M_in);
            axis_tmp = cell2mat(rho_num_cell(M_in ~= sym(0)));
            coef_tmp = M_in(M_in ~= sym(0));

            obj = PO();
            obj.axis = axis_tmp;
            obj.coef = coef_tmp;
            obj.spin_label = spin_label_cell;
            obj.bracket = zeros(size(coef_tmp));
            obj.basis = 'pol';

            obj = CombPO(obj);
            obj.logs = obj.txt;
        end
        % M2pol

        %% obj = M2pmz(M_in,spin_label_cell)
        function obj = M2pmz(M_in,spin_label_cell)
            if nargin < 2
                spin_label_cell = PO.spin_label_cell_default; % Default spin_label
            end
            obj = PO.M2pol(M_in,spin_label_cell);
            obj = pol2pmz(obj);
            obj.logs = obj.txt;
        end
        % M2pmz

        %% obj = M2xyz(M_in,spin_label_cell)
        function obj = M2xyz(M_in,spin_label_cell)
            if nargin < 2
                spin_label_cell = PO.spin_label_cell_default; % Default spin_label
            end
            obj = PO.M2pol(M_in,spin_label_cell);
            obj = pol2xyz(obj);
            obj.logs = obj.txt;
        end
        % M2xyz

    end % methods (Static)
end % classdef