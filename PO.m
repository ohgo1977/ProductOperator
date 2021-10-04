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
% classdef (InferiorClasses = {?sym}) PO 
% Overloading of some functions of sym class

% https://jp.mathworks.com/help/matlab/matlab_oop/custom-display-interface.html

    %%    
    properties (SetAccess = protected) % Read-Only from the Command Window
        axis        % Showing the status of axis direction for each spin.
                    % 0:E 1:x 2:y 3:z 4:p 5:m 6:a 7:b
                    % The column size corresponds to the number of spin types in the system.

        coef        % Coefficients of product operators (Symbolic).
                    % coef should not include the 2^(N-1) coefficient.

        spin_label  % Labels for spin1, 2, 3... stored in a cell. 
                    % Default: {'I' 'S' 'K' 'L' 'M'} defined in the method PO().

        basis       % String value to distinguish the basis-status in the calculations.
                    % 'xyz', 'pmz' or 'pol'
                    
    end

    %%
    properties (Access = protected) % No access from the Command Window
        bracket % Binary value to indicate cases with (a+b) or (a-b) type coefficient (1: yes, 0: no)

        SimplifySteps = 10
                % Number of steps used for simplify() in CombPO().
                % This prorperty can be changed from set_SimplifySteps().
    end

    %% 
    properties (Dependent) % Parameters depending on other parameters
        Ncoef       % The 2^(N-1) coefficient for N spin-1/2 system (Symbolic) 
        txt         % Text output of Product Operators (String)
        M           % Matrix Form
        coherence   % Populations (diagonal) and coherences (off-diagonal) of a density operator
    end
    
    %%
    properties (Constant = true) % Constant throughout the methods
        sqn = sym(1/2);% Spin Quantum Number
    end
    
    %%
    properties % No access limit
        disp = 1    % Control the display of the applied method on the monitor.
                    % 1: On, 2: Off
    end

    % Custom Display of properties. Need to learn more.
    % These methods worked but were slow.
    % Change classdef line.
    methods (Access = protected)
        % function header = getHeader(obj)
        %     if ~isscalar(obj)
        %         header = getHeader@matlab.mixin.CustomDisplay(obj);
        %     else
        %         % headerStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
        %         headerStr = getHeader@matlab.mixin.CustomDisplay(obj);
        %         M_txt = char(obj.M);
        %         id_tmp = strfind(M_txt,';');
        %         id_tmp = [0 id_tmp length(M_txt)];
        %         s_tmp = sprintf('M:');
        %         for ii = 1:length(id_tmp)-1
        %             s_tmp_tmp = M_txt(id_tmp(ii)+1:id_tmp(ii+1));
        %             s_tmp = sprintf('%s\n%s',s_tmp,s_tmp_tmp);
        %         end
        %         header = sprintf('%s\n%s',headerStr,s_tmp);
        %     end
        %  end
        function footer = getFooter(obj)
            if ~isscalar(obj)
                footer = getFooter@matlab.mixin.CustomDisplay(obj);
            else
                % footer = '';
                % for ii = 1:2
                %     if ii == 1
                %         M_txt = char(obj.M);
                %     elseif ii == 2
                %         M_txt = char(obj.coherence);
                %     end
                %     id_tmp = strfind(M_txt,';');
                %     id_tmp = [0 id_tmp length(M_txt)];
                %     if ii == 1
                %         footer = sprintf('%sM:',footer);
                %     else
                %         footer = sprintf('%s\n\ncoherence:',footer);
                %     end

                %     for jj = 1:length(id_tmp)-1
                %         s_tmp = M_txt(id_tmp(jj)+1:id_tmp(jj+1));
                %         footer = sprintf('%s\n%s',footer,s_tmp);
                %     end
                % end

                % https://www.mathworks.com/matlabcentral/answers/6940-save-disp-output-as-a-string-cell-variable
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
                propList.basis = obj.basis;
                propList.disp = obj.disp;
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
            elseif strcmp(obj.basis, 'pmz') || strcmp(obj.basis, 'pol')% Raising/Lowering  operator basis
                Ncoef_out = sym(ones(size(obj.axis,1),1));% Ncoef = 1
            end
        end % get.Ncoef
        
        %% txt_out = get.txt(obj)
        function txt_out = get.txt(obj)
            if isempty(find(obj.coef ~= sym(0),1))% If all coef values are zero.
                txt_out = '0';
            else
                txt_out = '';
                for ii = 1:length(obj.coef)
                    axis_tmp = obj.axis(ii,:);
                    Ncoef_tmp = obj.Ncoef(ii);
                    coef_tmp = obj.coef(ii);
                    bracket_tmp = obj.bracket(ii);

                    % Text of Product Operator
                    pt = axis2pt(obj,axis_tmp);
                    
                    % Remove '1' from single-type P.O. (N = 1 for 2^(N-1))
                    if ~strcmp(char(Ncoef_tmp),'1')
                        ptc = strcat(char(Ncoef_tmp),pt);
                    else
                        ptc = pt;
                    end

                    % Adjustment of sign and Creation of txt
                    subexpr = children(coef_tmp);
                    if sign(subexpr{end}) == -1 || sign(coef_tmp) == -1 || sign(subexpr{end}) == -1i || ... 
                    (length(subexpr) == 2 && sign(subexpr{1}) == -1 && sign(subexpr{2}) == 1i)
                        % Case of negative values: change the position of '-'.
                        % 1st condition: Symbols with negative sign such as -q, -1/2*q, etc..
                        % 2nd condition: symbolic negative values such as sym(-2).
                        % 3rd condition: symbols with negative imaginary such as -a*1i
                        % 4th condition: symbolic negative imaginary valuess such as -5*1i
                        if ~strcmp(char(coef_tmp),'-1')% if coef_tmp = sym(-1), it is not required to display -1 as a coeffcieint.
                            if bracket_tmp == 1
                                ptc = strcat(ptc,'*','(',char(-1*coef_tmp),')');% Add bracket for 'a-b'-type coefficient
                            else
                                ptc = strcat(ptc,'*',char(-1*coef_tmp));
                            end
                        end
                        txt_out = [txt_out,' ','-',' ',ptc];% Add Negative sign in the text

                    else % Symbols with positive sign such as q in addtion to symbolic positive numbers.
                        if ~strcmp(char(coef_tmp),'1')% if coef_tmp = sym(1), it is not required to display 1 as a coeffcieint.
                            if bracket_tmp == 1
                                ptc = strcat(ptc,'*','(',char(coef_tmp),')');% Add bracket for 'a+b'-type coefficient
                            else
                                ptc = strcat(ptc,'*',char(coef_tmp));
                            end            
                        end

                        if ~strcmp(char(coef_tmp),'0')% if coef_tmp ~= sym(0)
                            if ii == 1% In the case of 1st term, no need to add '+'.
                                txt_out = [txt_out,ptc];% No Positive sign for 1st term
                            else
                                txt_out = [txt_out,' ','+',' ',ptc];
                            end
                        end

                    end
                end
            end
        end  % get.txt

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
                % Mo is a symbolic class because coef, Ncoef are symbolic.
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
        end % get.M(obj)
        
        %% coherence_out = get.coherence(obj)
        function coherence_out = get.coherence(obj)   
            % Populations and coherences of a density operator
            spin_no = size(obj.axis,2);
            coherence_out = PO.rho_box(spin_no);
            coherence_out(obj.M == 0) = 0;
        end

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
                spin_label_cell = {'I' 'S' 'K' 'L' 'M'}; % Default spin_label
                end

                if length(spin_label_cell) < spin_no % Abort spin_label_cell is not large enough.
                error('the size of spin_label_cell must be same as or bigger than spin_no');
                end

                spin_label_cell = spin_label_cell(1:spin_no);% Adjust the size of spin_label_cell to spin_no. 

                for ii = 1:max(size(sp_cell))
                    sp = sp_cell{ii};

                    axis_tmp = zeros(1,spin_no);
                    % Initial state of axis_temp is 1/2E for all spins.
                    % if sp is not listed in spin_lable,
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
                    error('Error: Cartesian operator basis and Raising/Lowering operator basis should not be mixed!!');
                end

                obj = PO(); % spin_label is empty at this point
                obj.axis = axis_out;% 1:x, 2:y, 3:z, 4:p, 5:m, 0: no type assgined
                obj.coef = coef_out;% Coefficient for the product operator
                obj.spin_label = spin_label_cell;
                obj.bracket = bracket_out;% 1: put bracket if coefficient is a sum-form.
                obj.basis = basis_out;

                obj = CombPO(obj);
            end
        end % PO
        
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
        
            for ii = 1:length(IA)
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

            A_dummy = sym('A_dummy'); % Use "A" as a first charcter so that A_dummy comes before other alphabets.
            % syms A_dummy % Use "A" as a first charcter so that A_dummy comes before other alphabets.
                         % High time-cost; This line should be out of the loop. 
            dummy_p_mat = A_dummy*coef_out;
            char_A_dummy = char(A_dummy);

            for ii = 1:length(coef_out)
                coef_tmp = coef_out(ii);
                dummy_p = dummy_p_mat(ii); % dummy_p = A_dummy*coef_tmp
                char_dummy_p = char(dummy_p);
                id_tmp = strfind(char_dummy_p,char(A_dummy)) + length(char_A_dummy) + 1;
        
                if coef_tmp == sym(0)% Special case: coef_tmp = sym(0)
                    bracket_tmp = 0;
                else
                    if length(char_dummy_p) > id_tmp && strcmp(char_dummy_p(id_tmp),'(')
                        % 1st condition is the case with coef_tmp = sym(1)
                        % if dummy_p = -3*A_dummy*(...), 
                        % then the charcter at id_tmp of char_dummy_p should be '('.
                        bracket_tmp = 1;
                    else
                        bracket_tmp = 0;                   
                    end
                end
                bracket_out(ii) = bracket_tmp;
            end

            % Remove terms with 0 coefficients
            if isempty(find(coef_out ~= sym(0),1)) % Only zero values in coef_out, then set as 0*1/2E. 
                id_vec = 1;
                axis_out = zeros(1,size(axis_in,2));% Reset axis_out for 1/2E
                bracket_out = 0;
            else
                id_vec = find(coef_out ~= sym(0));
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

        end %CombPO

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
        end
        % pol2xyz

        %% obj = pmz2pol(obj)
        function obj = pmz2pol(obj)
            if ~strcmp(obj.basis,'pmz')
                error('The basis of the object should be pmz')
            end

            obj = pmz2xyz(obj);
            obj = xyz2pol(obj);
            % Future version: direct conversin from pmz to pol
        end
        % pmz2pol

        %% obj = pol2pmz(obj)
        function obj = pol2pmz(obj)
            if ~strcmp(obj.basis,'pol')
                error('The basis of the object should be pol')
            end

            obj = pol2xyz(obj);
            obj = xyz2pmz(obj);
            % Future version: direct conversin from pmz to pol
        end
        % pol2pmz

        %% obj = UrhoUinv(obj,H,q)
        function obj = UrhoUinv(obj,H,q)
            % obj = UrhoUinv(obj,H,q)
            % Calculation of the evolution of rho under q*H based on the cyclic
            % commutations. No matrix calculation is used.
            %
            % obj, H: PO class objects with xyz-basis.
            % The size of H.axis should be [1,n], i.e., only single term.
            % q: angle in radian (symbolic or double)
            
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

            if strcmp(obj.basis,'pmz') || strcmp(H.basis,'pmz') || strcmp(obj.basis,'pol') || strcmp(H.basis,'pol') 
                error('The basis of the object should be xyz')
            end

            mt = [ 0  3 -2; 
                  -3  0  1; 
                   2 -1  0];

            % Conversion of q from double to Symbolic
            q = sym(q);

            % Calculation of new density operator evolved under a Hamiltonian
            axis_new = [];
            coef_new = [];
            H_axis = H.axis;

            type_mask_mat = (obj.axis.*H_axis) ~= 0;% Check how many spin types get matched, matched: 1, unmatched: 0
            axis_diff_mat = obj.axis ~= H_axis;% Check the difference of the direction of each spin type
            axis_mask_mat = type_mask_mat.*axis_diff_mat;
            axis_mask_vec = sum(axis_mask_mat,2);

            for ii = 1:length(obj.coef)% For each term of rho
                rho_axis = obj.axis(ii,:);
                axis_mask = axis_mask_vec(ii);

                % type_mask_vec = (rho_axis.*H_axis) ~= 0;% Check how many spin types get matched, matched: 1, unmatched: 0
                % axis_diff_vec = rho_axis ~= H_axis;% Check the difference of the direction of each spin type
                % axis_mask_vec = type_mask_vec.*axis_diff_vec;
                % axis_mask = sum(axis_mask_vec);

                % The cases of axis_mask = 1 means
                % there is at least one match of spin types between H and rho
                % and
                % only one spin-type has a difference of axis directions between H
                % and rho.
                % The evolution of rho under H should be considered with
                % these conditions.
                
                axis_new = cat(1,axis_new,rho_axis);% Put original axis in new axis

                if  axis_mask == 1
                    axis_tmp = zeros(1,size(obj.axis,2));% This should be 0, i.e. Unit-operator.
                    sign_tmp = 1;

                    for jj = 1:length(rho_axis)
                        H_axis_tmp = H_axis(jj);
                        rho_axis_tmp = rho_axis(jj);

                        % determination of sine term
                        if H_axis_tmp ~= 0 && rho_axis_tmp ~= 0 % Use of Master Table                           
                           mtv_tmp = mt(H_axis_tmp,rho_axis_tmp); 

                        elseif H_axis_tmp ~= 0 && rho_axis_tmp == 0 % Ex. rho: Ix, H: 2IzSz => sine term: 2Iy(Sz)
                           mtv_tmp = H_axis_tmp;                    % (Sz) corresonds to H_axis_tmp
                               
                        elseif H_axis_tmp == 0 && rho_axis_tmp ~= 0 % Ex. rho: 2IzSz, H: Ix => sine term: 2Iy(Sz)
                           mtv_tmp = rho_axis_tmp;                  % (Sz) corresonds to rho_axis_tmp

                        else
                           mtv_tmp = 0; % No evolution. Does this case happen with axis_mask == 1?

                        end

                        if mtv_tmp ~= 0
                            axis_tmp(jj) = abs(mtv_tmp);
                            sign_tmp = sign_tmp*sign(mtv_tmp);
                        end
                    end
                    
                    sign_tmp = sign_tmp*sign(H.coef);% Direction of the rotation by the sign of H

                    if ~isempty(find(axis_tmp,1))% if axis_tmp is not [0 0 0]
                        axis_new = cat(1,axis_new,axis_tmp);
                        coef_new = cat(1,coef_new,obj.coef(ii)*[cos(q); sign_tmp*sin(q)]);
                    end

                else % axis_mask ~= 1, No evolution
                    coef_new = cat(1,coef_new,obj.coef(ii));
                end
            end

            obj.axis = axis_new;
            obj.coef = coef_new;
            obj.bracket = zeros(size(obj.coef));
            obj = CombPO(obj);
        end % UrhoUinv
        
        %% obj = pulse(obj,sp,ph,q)
        function [obj, id_sp] = pulse(obj,sp,ph,q)
            % obj = pulse(obj,sp,ph,q)
            % Calculation of the change of rho under a pulse.
            % obj: PO class object
            % sp: type of spin, character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % ph: phase of the pulse, character ('x','X','-x' etc.) or number (0,1,2,3)
            %     quadrature phase only.
            % q: flip angle in radian (symbolic or double)
            
            basis_org = 'xyz';
            if strcmp(obj.basis,'pmz')
                obj =  pmz2xyz(obj);% Conversion to xyz
                basis_org = 'pmz';
            elseif strcmp(obj.basis,'pol')
                obj =  pol2xyz(obj);% Conversion to xyz
                basis_org = 'pol';
            end

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
            H.bracket = 0;

            obj = UrhoUinv(obj,H,q);

            if strcmp(basis_org,'pmz')
                obj =  xyz2pmz(obj);
            elseif strcmp(basis_org,'pol')
                obj =  xyz2pol(obj);
            end

            if isa(q,'sym') == 1
                s_out = sprintf('Pulse: %s %s %s',spin_label_cell{id_sp},char(q),ph_tmp);
            else
                s_out = sprintf('Pulse: %s%4d%s',spin_label_cell{id_sp},round(q/pi*180),ph_tmp);
            end

            if obj.disp == 1
                fprintf(1,'%s\n',s_out);
                fprintf(1,'    %s\n',obj.txt);
            end
        end % pulse
        
        %% obj = simpulse(obj,sp_cell,ph_cell,q_cell)
        function obj = simpulse(obj,sp_cell,ph_cell,q_cell)
            % obj = simpulse(obj,sp_cell,ph_cell,q_cell)
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
            % rho_new = simpulse(rho,{'I' 'S'},{'x' 'y'},{pi/2 pi/2})
            
            obj_tmp = obj;
            spin_label_cell = obj.spin_label;
            for ii = 1:max(size(sp_cell))
                sp = sp_cell{ii};
                ph = ph_cell{ii};
                q = q_cell{ii};
             
                if contains(sp,'*')
                    if length(sp) == 1 % '*' applys pulse() to all spins.
                        id_vec = 1:max(size(spin_label_cell));
                    elseif length(sp) == 2 % 'I*' applys pulse() to I1, I2, ... if spin_label is set {'I1' 'I2' ...}
                        id_vec = find(contains(spin_label_cell,sp(1)));
                    end

                    for jj = id_vec
                        sp_tmp = spin_label_cell{jj};
                        obj_tmp = pulse(obj_tmp,sp_tmp,ph,q);%Run pulse each
                    end                   
                else % sp doesn't include '*'
                    obj_tmp = pulse(obj_tmp,sp,ph,q);%Run pulse each
                end
            end
            obj = obj_tmp;

        end % simpulse
        
        %% obj = cs(obj,sp,q)
        function [obj, id_sp] = cs(obj,sp,q)
            % obj = cs(obj,sp,q)
            % Calculation of the chemical shift evolution of rho.
            % obj: PO class object
            % sp: type of spin, character of spin ('I' or 'S' etc.) or
            % the order number of spin (1 for 'I', 2 for 'S' etc.).
            % q: flip angle (symbolic or double)

            basis_org = 'xyz';
            if strcmp(obj.basis,'pmz')
                obj =  pmz2xyz(obj);% Conversion to xyz
                basis_org = 'pmz';
            elseif strcmp(obj.basis,'pol')
                obj =  pol2xyz(obj);% Conversion to xyz
                basis_org = 'pol';
            end

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
            H.bracket = 0;% 1: put bracket if coefficient is a sum-form.

            obj = UrhoUinv(obj,H,q);

            if strcmp(basis_org,'pmz')
                obj =  xyz2pmz(obj);% Conversion to pmz
            end

            if isa(q,'sym') == 1
                s_out = sprintf('CS: %s %s',spin_label_cell{id_sp},char(q));
            else
                s_out = sprintf('CS: %s%4d',spin_label_cell{id_sp},round(q/pi*180));
            end

            if obj.disp == 1
                fprintf(1,'%s\n',s_out);
                fprintf(1,'    %s\n',obj.txt);
            end

        end % cs
        
        %% obj = simcs(obj,sp_cell,q_cell)
        function obj = simcs(obj,sp_cell,q_cell)
            % obj = simcs(obj,sp_cell,q_cell)

            obj_tmp = obj;
            spin_label_cell = obj.spin_label;
            for ii = 1:max(size(sp_cell))
                sp = sp_cell{ii};
                q = q_cell{ii};
             
                if contains(sp,'*')% Including 'I*' or '*'
                    if length(sp) == 1 % '*'
                        id_vec = 1:max(size(spin_label_cell));
                    elseif length(sp) == 2 % 'I*' 
                        id_vec = find(contains(spin_label_cell,sp(1)));
                    end

                    for jj = id_vec
                        sp_tmp = spin_label_cell{jj};
                        obj_tmp = cs(obj_tmp,sp_tmp,q);
                     end
                else % sp doesn't include '*'
                    obj_tmp = cs(obj_tmp,sp,q);
                end
             end
             obj = obj_tmp;
        end % simcs        
        
        %% obj = jc(obj,sp,q)
        function [obj, id_sp] = jc(obj,sp,q)
            % obj = jc(obj,sp,q)

            basis_org = 'xyz';
            if strcmp(obj.basis,'pmz')
                obj =  pmz2xyz(obj);% Conversion to xyz
                basis_org = 'pmz';
            elseif strcmp(obj.basis,'pol')
                obj =  pol2xyz(obj);% Conversion to xyz
                basis_org = 'pol';
            end

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
            H.bracket = 0;

            obj = UrhoUinv(obj,H,q);

            if strcmp(basis_org,'pmz')
                obj =  xyz2pmz(obj);% Conversion to pmz
            end

            if isa(q,'sym') == 1
                s_out = sprintf('JC: %s %s',sp_tmp,char(q));
            else
                s_out = sprintf('JC: %s%4d',sp_tmp,round(q/pi*180));
            end

            if obj.disp == 1
                fprintf(1,'%s\n',s_out);
                fprintf(1,'    %s\n',obj.txt);
            end

        end % jc 
       
        %% obj = simjc(obj,sp_cell,q_cell)
        function obj = simjc(obj,sp_cell,q_cell)
            % obj = simjc(obj,sp_cell,q_cell)
            obj_tmp = obj;
            for ii = 1:max(size(sp_cell))
               sp = sp_cell{ii};
               q = q_cell{ii};
               obj_tmp = jc(obj_tmp,sp,q);
            end
            obj = obj_tmp;
        end


        %% dispPOtxt(obj)
        function dispPOtxt(obj)
            % dispPOtxt(obj)
            % Display obj.txt
            fprintf(1,'    %s\n',obj.txt);
        end % dispPOtxt
       
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
                    if strcmp(char(obj.Ncoef(ii)),'1')
                        fprintf(1,'%4d%5s%-10s %s\n',ii,'',pt,char(obj.coef(ii)));
                    else
                        fprintf(1,'%4d%5s%-10s %s\n',ii,char(obj.Ncoef(ii)),pt,char(obj.coef(ii)));                   
                    end
                end
                fprintf(1,'\n');
        end % dispPO
       
        %% obj = pulse_phshift(obj,sp,ph,q)
        function obj = pulse_phshift(obj,sp,ph,q)
            % obj = pulse_phshift(obj,sp,ph,q)
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
            obj = cs(obj,sp,-ph);        % 1. Rotation -ph along Z axis
            obj = pulse(obj,sp,'x',q);   % 2. Rotation   q along X axis
            [obj, id_sp] = cs(obj,sp,ph);% 3. Rotation  ph along Z axis.
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

            if obj.disp == 1
                fprintf(1,'%s\n',s_out);
                fprintf(1,'    %s\n',obj.txt);
            end

        end % pulse_phshift
        
        %% obj = simpulse_phshift(obj,sp_cell,ph_cell,q_cell)
        function obj = simpulse_phshift(obj,sp_cell,ph_cell,q_cell)
            % obj = simpulse_phshift(obj,sp_cell,ph_cell,q_cell)

            obj_tmp = obj;
            spin_label_cell = obj.spin_label;
            for ii = 1:max(size(sp_cell))
                sp = sp_cell{ii};
                ph = ph_cell{ii};
                q = q_cell{ii};
             
                if contains(sp,'*')% Including 'I*' or '*'
                    if length(sp) == 1 % '*'
                        id_vec = 1:max(size(spin_label_cell));
                    elseif length(sp) == 2 % 'I*' 
                        id_vec = find(contains(spin_label_cell,sp(1)));
                    end

                    for jj = id_vec
                        sp_tmp = spin_label_cell{jj};
                        obj_tmp = pulse_phshift(obj_tmp,sp_tmp,ph,q);%Run pulse each
                     end                   
                else % sp doesn't include '*'
                    obj_tmp = pulse_phshift(obj_tmp,sp,ph,q);%Run pulse each
                end
            end
            obj = obj_tmp;

        end % simpulse_phshift

        %% obj = pfg(obj,G,gamma_cell)
        function obj = pfg(obj,G,gamma_cell)
            % obj = pfg(obj,G,gamma_cell)
            % applys pulse field gradient to all spins.
            % G is a strengh of the field and 
            % gamma_cell a cell array including gyromagnetic ratio of the spins.

            obj_tmp = obj;
            spin_label_cell = obj.spin_label;
            id_vec = 1:size(obj.axis,2);

            for jj = id_vec
                sp_tmp = spin_label_cell{jj};
                q = G*gamma_cell{jj};
                obj_tmp = cs(obj_tmp,sp_tmp,q);
            end
            obj = obj_tmp;
        end % pfg        
        
        %% a0_V = SigAmp1(obj,sp,phR)
        function [a0_V,rho_V] = SigAmp1(obj,sp,phR)
            % a0_V = SigAmp1(obj,sp,phR)
            % Calculation of initial signal amplitudes (t=0) in the equation
            % s(t) = 2*i*(rho[-b](t) + rho[-a](t) + rho[b-](t) + rho[a-](t))*exp(-i*phrec)
            % Spin Dynamics (2nd Ed.), p.379.
            %
            % Related topics: Spin Dynamics (2nd Ed.), p.262, p. 287, p. 371, p.379, pp.608-610.
            %
            % Example
            % a0_V = SigAmp(rho,'S','y')
            % a0_V = SigAmp(rho,'IS',0)
            % a0_V = SigAmp(rho,[1 2],0)

                spin_label_cell = obj.spin_label;

                if isa(sp,'double')
                    for ii = 1:length(sp)
                        sp_tmp = sp(ii); % double                   
                        sp_tmp = spin_label_cell{sp_tmp};% char 
                        sp_m = [sp_tmp 'm'];% sp_m = 'ImSm ...'.

                        ObsPO = PO(size(obj.axis,2),{sp_m},{sym(1)},spin_label_cell);
                        if ii == 1
                            obsPO_M = ObsPO.M;% Create obsPO_M
                        else
                            obsPO_M = obsPO_M + ObsPO.M;
                        end
                    end
                elseif isa(sp,'char')
                    ii_int = 0;
                    for ii = 1:length(spin_label_cell)
                        spin_label_tmp = spin_label_cell{ii};
                        if contains(sp,spin_label_tmp)
                            ii_int = ii_int + 1;
                            sp_tmp = spin_label_tmp;% char 
                            sp_m = [sp_tmp 'm'];% Im, Sm, ... .
        
                            ObsPO = PO(size(obj.axis,2),{sp_m},{sym(1)},spin_label_cell);
                            if ii_int == 1
                                obsPO_M = ObsPO.M;% Create obsPO_M
                            else
                                obsPO_M = obsPO_M + ObsPO.M;
                            end
                        end
                    end
                end

                a0_M = obj.M.*obsPO_M;            
                % This should be Hadamard product
                % rho.M.*(Im.M + Sm.M + ...) extracts only (-1)-quantum coherence components in rho, 
                % i.e., (Im.M + Sm.M + ...) works as a mask.
                
                obsPO_V = reshape(obsPO_M,1,numel(obsPO_M));
                id_tmp = obsPO_V ~= sym(0);

                a0_V = reshape(a0_M,1,numel(a0_M));
                % id_tmp = a0_V ~= sym(0);
                a0_V = a0_V(id_tmp);
                a0_V = 2*1i*PO.rec_coef(phR)*a0_V;
                
                rho = obj.coherence;
                rho_V = reshape(rho,1,numel(rho));
                rho_V = rho_V(id_tmp);

                if obj.disp == 1
                    ph_s = PO.ph_num2str(phR);
                    fprintf(1,'phRec: %2s\n',ph_s);
                end                
        end % SigAmp1

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

                ObsPO = PO(size(obj.axis,2),{sp_m},{sym(1)},spin_label_cell);
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

            a0_V = 2*1i*PO.rec_coef(phR)*a0_V;

            rho = obj.coherence;

            rho_V = reshape(rho,1,numel(rho));
            rho_V = rho_V(id_tmp);

            if obj.disp == 1
                ph_s = PO.ph_num2str(phR);
                fprintf(1,'phRec: %2s\n',ph_s);
            end

        end % SigAmp
       
       %% pt = axis2pt(obj,axis_tmp)
       function pt = axis2pt(obj,axis_tmp)
           % pt = axis2pt(obj,axis_tmp)
           % axis_tmp is a row vector from rho.axis.
           % obj is necessary to get spin_label from it.
           %
           % Example:
           % axis_tmp = [1 1], pt = 'IxSx'. Note that pt doesn't include '2' for '2IxSx'.
            pt = '';
            if isempty(find(axis_tmp, 1))
                pt = 'E';
            else
                for jj = 1:length(axis_tmp)
                    axis_v = axis_tmp(jj);
                    st = obj.spin_label{jj};

                    if axis_v ~=0
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
                        pt = strcat(pt,st,at);
                    end            
                end
            end
       end % axis2pt

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

            id_in = obj.findterm(id_in);
            obj.axis(id_in,:) = [];
            obj.coef(id_in,:) = [];
            obj = CombPO(obj);
        end % delPO

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

            id_in = obj.findterm(id_in);
            obj.axis = obj.axis(id_in,:);
            obj.coef = obj.coef(id_in,:);

            obj = CombPO(obj);
        end % selPO

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
        end % findcoef

        %% obj3 = commutator(obj1, obj2)
        function obj3 = commutator(obj1, obj2)
            % obj3 = commutator(obj1, obj2)
            % Commutation between obj1 and obj2.
            % If obj1 and obj2 don't commute, obj3 is 0
            % if [A,B] = iC, then B ==> B*cos(q) + C*sin(q) under A. 
            obj3 = obj1*obj2 - obj2*obj1;
        end % commute

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

       function obj = changeProp(obj, PropertyName, para_in)
            eval(['obj.',PropertyName, '= para_in;']);
       end
       % changeProp

        %% obj = set_SimplifySteps(ojb,new_v)
        function obj = set_SimplifySteps(obj,new_v)
            % obj = set_SimplifySteps(ojb,new_v)
            % Change the property SimplifySteps to a new value.
            obj.SimplifySteps = new_v;
        end
        % set_SimplifySteps

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
                    if strcmp(obj1.basis,'xyz')
                        obj1 = xyz2pmz(obj1);
                    end

                    if strcmp(obj2.basis,'xyz')
                        obj2 = xyz2pmz(obj2);
                    end
                    branch_id = 1;

                % xyz + pol => pol
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz'))
                    if strcmp(obj1.basis,'xyz')
                        obj1 = xyz2pol(obj1);
                    end

                    if strcmp(obj2.basis,'xyz')
                        obj2 = xyz2pol(obj2);
                    end
                    branch_id = 2;

                % pmz + pol => pol
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    if strcmp(obj1.basis,'pmz')
                        obj1 = pmz2pol(obj1);
                    end

                    if strcmp(obj2.basis,'pmz')
                        obj2 = pmz2pol(obj2);
                    end                    
                    branch_id = 3;
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
                    if strcmp(obj1.basis,'xyz')
                        obj1 = xyz2pmz(obj1);
                    end

                    if strcmp(obj2.basis,'xyz')
                        obj2 = xyz2pmz(obj2);
                    end
                    branch_id = 1;

                % xyz - pol => pol
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz'))
                    if strcmp(obj1.basis,'xyz')
                        obj1 = xyz2pol(obj1);
                    end

                    if strcmp(obj2.basis,'xyz')
                        obj2 = xyz2pol(obj2);
                    end
                    branch_id = 2;

                % pmz - pol => pol
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    if strcmp(obj1.basis,'pmz')
                        obj1 = pmz2pol(obj1);
                    end

                    if strcmp(obj2.basis,'pmz')
                        obj2 = pmz2pol(obj2);
                    end                    
                    branch_id = 3;
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
        end
        % minus

        %% obj = uminus(obj) 
        function obj = uminus(obj)
            obj.coef = -1*obj.coef;
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

                % xyz * pmz
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pmz')) || (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'xyz'))
                    if strcmp(obj1.basis,'pmz')
                        obj1 = pmz2xyz(obj1);
                    end

                    if strcmp(obj2.basis,'pmz')
                        obj2 = pmz2xyz(obj2);
                    end
                    basis_org = 'pmz';
                    branch_id = 2;

                % xyz * pol
                elseif (strcmp(obj1.basis,'xyz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'xyz'))
                    if strcmp(obj1.basis,'pol')
                        obj1 = pol2xyz(obj1);
                    end

                    if strcmp(obj2.basis,'pol')
                        obj2 = pol2xyz(obj2);
                    end
                    basis_org = 'pol';
                    branch_id = 3;

                % pmz * pmz
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pmz'))
                    obj1 = pmz2xyz(obj1);
                    obj2 = pmz2xyz(obj2);
                    basis_org = 'pmz';
                    branch_id = 4;

                % pmz * pol
                elseif (strcmp(obj1.basis,'pmz') && strcmp(obj2.basis,'pol')) || (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pmz'))
                    if strcmp(obj1.basis,'pol')
                        obj1 = pol2xyz(obj1);
                    elseif strcmp(obj1.basis,'pmz')
                        obj1 = pmz2xyz(obj1);
                    end

                    if strcmp(obj2.basis,'pol')
                        obj2 = pol2xyz(obj2);
                    elseif strcmp(obj2.basis,'pmz')
                        obj2 = pmz2xyz(obj2);
                    end
                    basis_org = 'pol';
                    branch_id = 5;

                % pol * pol
                elseif (strcmp(obj1.basis,'pol') && strcmp(obj2.basis,'pol'))
                    obj1 = pol2xyz(obj1);
                    obj2 = pol2xyz(obj2);
                    basis_org = 'pol';
                    branch_id = 6;

                end

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

                % Conversion from xyz to pmz
                if strcmp(basis_org, 'pmz')
                    obj = xyz2pmz(obj);
                elseif strcmp(basis_org, 'pol')
                    obj = xyz2pol(obj);
                end
            else
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
            end
        end
        % mtimes        

        %% obj = mrdivide(obj1, obj2)
        function obj = mrdivide(obj1, obj2) %obj1/obj2
            if isa(obj2,'double')||isa(obj2,'sym') ||isa(obj2,'char') % obj1/a
                obj_base = obj1;
                coef_tmp = sym(obj2);
            else          
                error('PO-class object cannot the be the divisor!')
            end
            coef_new = obj_base.coef/coef_tmp;
            obj_base.coef = coef_new;
            obj = CombPO(obj_base);
        end
        % mrdivide        
    end % methods
    
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
        end % phmod
       
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
        end % ph_num2st
       
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
        end % rec_coef              
       
        %% M_out = axis2M(axis_v, sqn)
        function M_out = axis2M(axis_v, sqn)
            % M_out = axis2M(axis_v, sqn)
            % axis_v is a value for PO.axis, i.e., 0, 1, 2, 3, 4, 5, 6 and 7.
            % sqn is spin quantum number (Symbolic)
            % Currently, only sqn = sym(1/2) works.
            % M_out is 1/2*Pauli matrix (Symbolic)
            % (axis_v = 0 for s0 (E), 1 for s1(sx), 2 for s2(sy) and 3 for s3(sz))

            if nargin == 1
                sqn = sym(1/2);
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
            %   M_out = sym(M_tmp);
            M_out = M_tmp; % Output as double. It is converted to sym in get.M().
        end % axis2M

        %% rho = rho_box(n)
        function rho = rho_box(n)
            % Calculation of "box notation" of the density oeprator, rho.
            % a: alpha state
            % b: beta state
            % m: a => b coherence
            % p: b => a coherence
            % Spin Dynamics p. 470, p. 260, p. 160
            n_s = 2^n;
            dec = n_s-1:-1:0;
        
            bin_mat = de2bi(dec,n,'left-msb');
            % bin_mat = dec2bin(dec,n) - '0';
            % % a: 1, b: 0
            % % a a a = 1 1 1
            % % a a b = 1 1 0
            % % a b a = 1 0 1
            % % a b b = 1 0 0
            % % b a a = 0 1 1
            % % b a b = 0 1 0
            % % b b a = 0 0 1
            % % b b b = 0 0 0
        
            rho_cell = cell(n_s,n_s);
            % rho(c,r) = <c|rho|r>
            % Since the spin state should be read from right to left,
            % rho(c,r) correspond to |r> => |c>
            % Spin Dynamics, p. 470.

            for ii = 1:n_s% row
                for jj = 1:n_s% column
                    r_vec = bin_mat(jj,:);%|r>
                    c_vec = bin_mat(ii,:);%|c>
                    rho_vec = c_vec - r_vec;% To know |r> => |c>  
                    % rho_vec take -1, 0, 1.
                    % -1: m and 1:p
                    % 0: no change from |c> to |r>  
                    
                    rho_tmp = char(double('b')*ones(1,n));% 'bbbb...'
                    rho_tmp(r_vec == 1) = 'a';            % 'abaa...' correspond to r_vec
                    rho_tmp(rho_vec == 1) = 'p';          % 'pbaa...'
                    rho_tmp(rho_vec == -1) = 'm';         % 'pbam...'

                rho_cell{ii,jj} = rho_tmp;
                end
            end
            rho = sym(rho_cell);

        end
        % rho = rho_box(n)

        %% create(spin_label_cell)
        function create(spin_label_cell, symcoef_switch)
            %% Spin Operators
            spin_no = length(spin_label_cell);
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
                    assignin('base',sp,obj);
                end
            end
            obj = PO(spin_no,{'1'},{1},spin_label_cell);
            assignin('base','hE',obj);

            if nargin == 1
                symcoef_switch = 'on';
            end

            switch symcoef_switch
                case {'y','on'}
                    PO.symcoef(spin_label_cell)
            end
        end
        % create

        %% symcoef(spin_label_cell)
        function symcoef(spin_label_cell,add_cell)
            % Create pre-set Symbolic Coefficients

            spin_no = length(spin_label_cell);
            % frequency o
            for ii = 1:spin_no
                sp = spin_label_cell{ii};

                varname = ['o' sp];
                assignin('base',varname,sym(varname));

                varname = ['o' sp(1)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end

                varname = ['o' sp(end)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end

                varname = ['o' num2str(ii)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end
            end

            % coupling J
            for ii = 1:spin_no
                for jj = 1:spin_no
                    if ii < jj
                        sp_ii = spin_label_cell{ii};
                        sp_jj = spin_label_cell{jj};

                        varname = ['J' sp_ii sp_jj];
                        assignin('base',varname,sym(varname));

                        varname = ['J' sp_ii(1) sp_jj(1)];
                        if exist(varname,'var') == 0
                            assignin('base',varname,sym(varname));
                        end

                        varname = ['J' sp_ii(end) sp_jj(end)];
                        if exist(varname,'var') == 0
                            assignin('base',varname,sym(varname));
                        end

                        varname = ['J' num2str(ii) num2str(jj)];
                        if exist(varname,'var') == 0
                            assignin('base',varname,sym(varname));
                        end
                    end
                end
            end

            % gyromagnetic ratio g
            for ii = 1:spin_no
                sp = spin_label_cell{ii};

                varname = ['g' sp];
                assignin('base',varname,sym(varname));

                varname = ['g' sp(1)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end

                varname = ['g' sp(end)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end

                varname = ['g' num2str(ii)];
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end
            end

            symcoef_list = {'a' 'b' 't1' 't' 'q' 'B' 'd' 'G'};
            if nargin == 1
                add_cell = {};
            end
            symcoef_list = [symcoef_list add_cell];
            for ii = 1:length(symcoef_list)
                varname = symcoef_list{ii};
                if exist(varname,'var') == 0
                    assignin('base',varname,sym(varname));
                end
            end
        end
        % symcoef
      
    end % methods (Static)
end % classdef