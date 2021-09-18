% ------------------------------------------------------------------------
% Class       : PO
% Description : Functions for Product Operator Formalism of spin-1/2
% Requirement : MATLAB Symbolic Math Toolbox
% Developer   : Dr. Kosuke Ohgo
% ------------------------------------------------------------------------
%
% Please read the manual (PO_Manual.docx or PO_Manual.pdf) for details.
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
classdef PO    
    properties
        axis        % Showing the status of axis direction for each spin.
                    % 1:x 2:y 3:z 4:p 5:m 0:1/2E
                    % The size of column corresponds to the number of spins.

        coef        % Coefficients of terms (Symbolic).
                    % coef doesn't include the 2^(N-1) coefficient.

        spin_label  % Labels for spin1, 2, 3... stored in a cell. Default: {'I' 'S' 'K' 'L' 'M'}

        disp = 1    % Control the display of the applied method on the monitor. 1: On, 2: Off

    end
    
    properties (Access = protected)
        bracket % Binary value to indicate cases with (a+b) or (a-b) type coefficient (1: yes, 0: no)

        basis   % String value to distinguish the basis-status in the calculations ('xyz' or 'pmz')
                % Ideally it should be visible to users (but still protected).
                % dispProp(obj,'basis') can be used to check this property.

        simplifystep = 10
                % Number of steps used for simplify() in CombPO().
                % This prorperty can be changed from set_simplifystep().
    end
    
    properties (Dependent)
        % Learn more about Dependent and getter function.
        % https://www.mathworks.com/help/matlab/matlab_oop/property-get-methods.html

        Ncoef       % The 2^(N-1) coefficient for N spin-1/2 system (Symbolic) 
        txt         % Text output of Product Operators (String)
        M           % Matrix Form
        coherence   % Populations (diagonal) and coherences (off-diagonal) of a density operator
    end
    
    properties (Constant = true)
        sqn = sym(1/2);% Spin Quantum Number
     %   url = 'https://github.com/ohgo1977/ProductOperator';
     %   version = '1.0.0';% Major.Minor.Patch

    end
        
    methods
        %% Ncoef_cnst = get.Ncoef(obj)
        function Ncoef_cnst = get.Ncoef(obj)
            if strcmp(obj.basis, 'xyz')% Cartesian operator basis
                Ncoef_cnst = sym(2.^(sum((obj.axis~=0),2)-1));
                % 2^(Ns-1) 
                % Examples 
                % IxSx   => Ns = 2 => Ncoef = 2
                % IxSxKx => Ns = 3 => Ncoef = 4
            elseif strcmp(obj.basis, 'pmz')% Raising/Lowering  operator basis
                Ncoef_cnst = sym(ones(size(obj.axis,1),1));
            end
        end % get.Ncoef
        
        %% txt_out = get.txt(obj)
        function txt_out = get.txt(obj)
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
        end  % get.txt

        %% M_out = get.M(obj)
        function M_out = get.M(obj)
            % Create a Matrix representation.
            for ii = 1:size(obj.axis,1)
                axis_tmp = obj.axis(ii,:);% should include the case with 0: 1/2E

                for jj = 1:length(axis_tmp)
                    Ma = obj.axis2M(axis_tmp(jj),obj.sqn); 
                    if jj == 1
                        M_tmp = Ma;
                    else
                        M_tmp = kron(M_tmp, Ma);% Kronecker tensor product
                    end
                end
                
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

        %% obj = PO(spin_no,sp_cell,symcoef_cell,spin_label_cell)
        function obj = PO(spin_no,sp_cell,symcoef_cell,spin_label_cell) 
            % obj = PO(spin_no,sp_cell,symcoef_cell,spin_label_cell)
            % This is a class constructor, obj should not be involved in argument.
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
            % symcoef_cell: Assignment of coefficients for sp_cell.
            % If only spin_no and sp_cell are input, sym(1) is automatically assigned for each term.
            % Example
            % sym q
            % {cos(q) sin(q)}
            % Do not include 2^(N-1) coefficients in symcoef_cell.
            %
            % spin_label_cell: Labels for spin types
            % If spin_label_cell is not input, {'I' 'S' 'K' 'L' 'M'} is used.
            % The size of spin_label_cell should be equal or larger than spin_no.
            % Example
            % spin_label_cell = {'I1' 'I2' 'I3' 'S4'}
            % In this case, sp_cell should be
            % sp_cell = {'I1x' 'I1yS4z'}
            % Note: there should be no overlap of strings among the members of spin_label_cell.
            % For example, Both 'I1' and 'I12' have 'I1'. This causes an error.
            %

            if nargin > 0 % This condition allows to create an empty PO object by obj = PO();
                axis_out = [];
                coef_out = [];
                bracket_out = [];

                 if nargin <= 3
                    spin_label_cell = {'I' 'S' 'K' 'L' 'M'}; % Default spin_label
                  end

                  if length(spin_label_cell) < spin_no % Abort spin_label_cell is not large enough.
                    error('the size of spin_label_cell must be same as spin_no');
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
                                otherwise, phase_id = 0;
                            end
                            axis_tmp(id_tmp) = phase_id;
                        end
                    end

                   axis_out = cat(1,axis_out,axis_tmp);

                   if nargin > 2
                       symcoef = symcoef_cell{ii};
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


                if isempty(find(axis_out == 4,1)) && isempty(find(axis_out == 5,1))% Cartesian operator basis
                    basis_out = 'xyz';
                elseif isempty(find(axis_out == 1,1)) && isempty(find(axis_out == 2,1))% Raising/Lowering operator basis
                    basis_out = 'pmz';
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
            % Combine coeffcieints of same type of terms in a PO-class object.
            % Also detect coefficients in which parentheses should be added and put a flag (obj.bracket).
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

            coef_out = [];
            axis_out = [];
            bracket_out = [];

            for ii = 1:length(IA)
               IA_tmp = IA(ii);
               IC_tmp = find(IC == IC(IA_tmp));% IMPORTANT!!
                % IA_tmp is an ID of a unique row in axis_in.
                % Then find(IC == IC(IA_tmp)) searches IDs of rows in axis_in corresponding to the unique row.
                % Example
                % IA_tmp = IA(1) = 2. 
                % IC(IA_tmp) = IC(2) = 1. 
                % Then find(IC == IC(IA_tmp)) = find(IC == 1) = [2 5];
                 
               axis_tmp = axis_in(IA_tmp,:);

               coef_tmp = sum(obj.coef(IC_tmp));
               coef_tmp = simplify(coef_tmp, 'Steps', obj.simplifystep);
      
               syms dummy_c
               dummy_p = dummy_c*coef_tmp;
               
               char_dummy_p = char(dummy_p);
               char_dummy_c = char(dummy_c);
               id_tmp = strfind(char_dummy_p,char(dummy_c))+length(char_dummy_c)+1;

               if coef_tmp == sym(0)% Special case: coef_tmp = sym(0)
                   bracket_tmp = 0;
               else
                   if length(char_dummy_p) > id_tmp && strcmp(char_dummy_p(id_tmp),'(')
                       % if dummy_p is -a*dummy_c*(...), 
                       % then the charcter at id_tmp of char_dummy_p should be '('.
                       bracket_tmp = 1;
                   else
                       bracket_tmp = 0;                   
                   end
               end

               axis_out = cat(1,axis_out,axis_tmp);
               coef_out = cat(1,coef_out,coef_tmp);
               bracket_out = cat(1,bracket_out,bracket_tmp);
            end

            % Remove terms with 0 coefficients
            id_vec = [];
            for ii = 1:length(coef_out)
                if ~strcmp(char(coef_out(ii)),'0')
                    id_vec = cat(1,id_vec,ii);
                end
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
        function obj = xyz2pmz(obj);
            % obj = xyz2pmz(obj)
            % conversion from Cartesian operator basis to lowring/raising operator basis

            if strcmp(obj.basis,'pmz')
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
                if isempty(find(axis_tmp == 1,1)) && isempty(find(axis_tmp == 2,1))% Iz, 2IzSz, 4IzSzKz,...
                    axis_out_tmp = axis_tmp;
                    coef_out_tmp = 1;
                else
                    % Conversion from Ix and Iy to Ip and Im
                    xn = length(find(axis_tmp == 1)); % Example: IxSyKxMz =>(Ip + Im)(Sp - Sm)(Kp + Km)Mz
                    yn = length(find(axis_tmp == 2)); % xn =2, yn = 1 => 2^(2+1) = 8 terms.
                    xyn = 2^(xn + yn);
                    axis_out_tmp = repmat(axis_tmp,xyn,1); % repmat([1 2 1 3],8,1)
                    axis_out_tmp(axis_out_tmp == 1) = 0;   % Remove x operators
                    axis_out_tmp(axis_out_tmp == 2) = 0;   % Remove y operators
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

                            c_tmp = zeros(size(bin_vec));
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
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);
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
            obj2.bracket = bracket_out;% change the property of bracket later.
            obj2.basis = 'pmz';
            obj = CombPO(obj2);
        end
        % xyz2pmz

        % obj = pmz2xyz(obj)
        function obj = pmz2xyz(obj);
            % obj = xyz2pmz(obj)
            % conversion from lowring/raising operator basis to Cartesian operator basis.

            if strcmp(obj.basis,'xyz')
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
                    xn = length(find(axis_tmp == 4));
                    yn = length(find(axis_tmp == 5));
                    zn = length(find(axis_tmp == 3));
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

                            c_tmp = zeros(size(bin_vec));
                            if axis_v == 4
                                c_tmp(bin_vec == 1) = 1;
                                c_tmp(bin_vec == 2) = 1i;
                            elseif axis_v == 5
                                c_tmp(bin_vec == 1) = 1;
                                c_tmp(bin_vec == 2) = -1i;
                            end
                            coef_out_tmp = coef_out_tmp.*c_tmp;
                        end
                    end
                end
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii)*(1/2)^(xn + yn +zn - 1);
                % From 1*IpSpIz*1, IxSxIz IxSyIz etc. are created.
                % Then Ncoef for IxSxIz in xyz-basis is calculated as 4 automatically. 
                % To compensate 4 in Ncoef, it is necessary to apply (1/2)*(xn + yn + zn -1) to a new coef.
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
            obj2.bracket = bracket_out;% change the property of bracket later.
            obj2.basis = 'xyz';
            obj = CombPO(obj2);
        end
        % pmz2xyz

        %% obj = UrhoUinv(obj,H,q)
        function obj = UrhoUinv(obj,H,q)
            % obj = UrhoUinv(obj,H,q)
            % Calculation of the evolution of rho under H based on the cyclic
            % commutations. No matrix calculation is used.
            % obj, H: PO class objects
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

            if strcmp(obj.basis,'pmz') || strcmp(H.basis,'pmz')
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
            for ii = 1:length(obj.coef)% For each term of rho
                rho_axis = obj.axis(ii,:);
                H_axis = H.axis;
                type_mask_vec = (rho_axis.*H_axis)~=0;% Check how many spin types get matched, matched: 1, unmatched: 0
                axis_diff_vec = rho_axis ~= H_axis;% Check the difference of the direction of each spin type
                axis_mask_vec = type_mask_vec.*axis_diff_vec;
                axis_mask = sum(axis_mask_vec);
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

                else
                    coef_new = cat(1,coef_new,obj.coef(ii));% No evolution.
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
            H.axis = axis_tmp;% 1:x, 2:y, 3:z, 0: no type assgined
            H.coef = coef_tmp;% Coefficient for the product operator
            H.bracket = 0;% 1: put bracket if coefficient is a sum-form.

            obj = UrhoUinv(obj,H,q);

            if strcmp(basis_org,'pmz')
                obj =  xyz2pmz(obj);% Conversion to pmz
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
        function obj = jc(obj,sp,q)
            % obj = jc(obj,sp,q)

            basis_org = 'xyz';
            if strcmp(obj.basis,'pmz')
                obj =  pmz2xyz(obj);% Conversion to xyz
                basis_org = 'pmz';
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
            elseif isa(sp,'char')
                for ii = 1:length(spin_label_cell)
                    spin_label_tmp = spin_label_cell{ii};
                    if contains(sp,spin_label_tmp)
                        id_tmp = ii;
                        axis_tmp(id_tmp) = 3;%Z direction
                    end
                end
                sp_tmp = sp;
            end
            
            H = PO();
            H.axis = axis_tmp;% 1:x, 2:y, 3:z, 0: no type assgined
            H.coef = sym(1);% Coefficient for the product operator
            H.bracket = 0;% 1: put bracket if coefficient is a sum-form.

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
            % expm(-1i*q*(Ix*cos(ph) + Iy*sin(ph)))*rho*expm(1i*q*(Ix*cos(ph) + Iy*sin(ph))) 
            % is equivalent to 
            % expm(-1i*ph*Iz)*expm(-1i*q*Ix)*expm(-1i*-ph*Iz)*rho*expm(1i*ph*Iz)*expm(1i*q*Ix)*expm(1i*-ph*Iz)
            % However,
            % expm(-1i*q*(Ix*cos(ph) + Iy*sin(ph)))*rho*expm(1i*q*(Ix*cos(ph) + Iy*sin(ph))) 
            % is not equivalent to
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
                    s_out = sprintf('Pulse: %s%4d%s',spin_label_cell{id_sp},round(q/pi*180),char(ph));
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
        
       %% a0_V = SigAmp(obj,sp,phR)
       function [a0_V,rho_V] = SigAmp(obj,sp,phR)
           % a0_V = SigAmp(obj,sp,phR)
           % Calculation of initial signal amplitudes in
           % the equation
           % s(t) = 2*i*(rho[-b](t) + rho[-a](t) + rho[b-](t) + rho[a-](t))*exp(-i*phrec)
           % Spin Dynamics (2nd Ed.), p.379.
           % Related topics: Spin Dynamics (2nd Ed.), p.262, p. 287, p. 371, p.379, pp.608-610.
           %
           % Example
           % a0_V = SigAmp(rho,'S','y')
           % a0_V = SigAmp(rho,'IS',0)
           % a0_V = SigAmp(rho,[1 2],0)

            spin_label_cell = obj.spin_label;

            % Observe Mx by Ix + Sx + ...
            if isa(sp,'double')
                for ii = 1:length(sp)
                    sp_tmp = sp(ii); % double                   
                    sp_tmp = spin_label_cell{sp_tmp};% char 
                    sp_x = [sp_tmp 'x'];% Ix, Sx, ... .

                    ObsPO = PO(size(obj.axis,2),{sp_x});
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
                        sp_x = [sp_tmp 'x'];% Ix, Sx, I1x, I2x,... .
    
                        ObsPO = PO(size(obj.axis,2),{sp_x});
                        if ii_int == 1
                            obsPO_M = ObsPO.M;% Create obsPO_M
                        else
                            obsPO_M = obsPO_M + ObsPO.M;
                        end
                    end
                end
            end
            
            a0_M = obj.M.*(2*tril(obsPO_M));
            % This should be Hadamard product
            % Do not describe as a0_M = obj.M.*2*tril(obsPO_M)!!
            % This description returns 2*obj.M*tril8obsPO_M.
            %
            % tril(Ix.M + Sx.M + ...) is I-.M + S-.M + ... where I-, S-, ... are down-shift operators.
            % Positions of the non-zero components (this case 1) in down-shift oeperators correspond to that of (-1)-quantum coherences in rho.
            % Thus, rho.M.*(2*tril(Ix.M + Sx.M + ...)) extracts only (-1)-quantum coherence components in rho, 
            % i.e., tril(Ix.M + Sx.M + ...) works as a mask.
            %
            % Since pmz basis operators can be used, this code may be written.
            
            a0_V = reshape(a0_M,1,numel(a0_M));
            id_tmp = a0_V ~= sym(0);
            a0_V = a0_V(id_tmp);
            a0_V = 2*1i*PO.rec_coef(phR)*a0_V;
            
            % syms rho [size(a0_M)] => rho1_2, rho2_2, ....
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
                                case '*', phase_id = [1;2;3;4;5];% Wildcard
                                otherwise, phase_id = 0;
                            end

                            if length(phase_id) == 1
                                    axis_tmp(:,id_tmp) = phase_id;
                            elseif length(phase_id) == 5% Wildcard for phase
                                rate_tmp = size(axis_tmp,1);%
                                axis_tmp = repmat(axis_tmp,5,1);% Expand the row size of axis_tmp 5 times. 
                                phase_id_vec = [1*ones(rate_tmp,1);...
                                                2*ones(rate_tmp,1);...
                                                3*ones(rate_tmp,1);...
                                                4*ones(rate_tmp,1);...
                                                5*ones(rate_tmp,1)];
                                % phase_id_vec = [1 1 1... 2 2 2... 3 3 3...4 4 4...5 5 5...]'
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
            obj3 = UrhoUinv(obj2,obj1,pi/2);

            if isempty(find(obj1.M - obj3.M ~= 0, 1))% obj1 == obj3
                % obj3 = PO(size(obj1.axis,2),{'1'},{sym(0)});
                fprintf(1,'They commutate\n')
            else
%                 obj3_tmp.coef = sym(1i);
                obj3.coef = obj3.coef*1i;
                obj3 = CombPO(obj3);
            end

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

    %    function obj = changeProp(obj, PropertyName, para_in)
    %         eval(['obj.',PropertyName, '= para_in;']);
    %    end
    %    % changeProp

        %% obj = set_simplifystep(ojb,new_v)
        function obj = set_simplifystep(obj,new_v)
            % obj = set_simplifystep(ojb,new_v)
            % Change the property simplifystep to a new value.
            obj.simplifystep = new_v;
        end
        % set_simplifystep(ojb,new_v)

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
        % the coef = exp(-1i*ph) = cos(ph) -1i*sin(ph)
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
        % axis_v is a value for PO.axis, i.e., 0, 1, 2, 3
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
              end
          end
          M_out = sym(M_tmp);
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
     
        rho = sym(zeros(n_s,n_s));
        % rho(c,r) = <c|rho|r>
        % Since the spin state should be read from right to left,
        % rho(c,r) correspond to |c> => |r>
        % Spin Dynamics, p. 470.
        for ii = 1:n_s% row
            for jj = 1:n_s% column
                r_vec = bin_mat(ii,:);%|r>
                c_vec = bin_mat(jj,:);%|c>
                rho_vec = r_vec - c_vec;% To know |c> => |r>  
                
                rho_tmp = '';
                for kk = 1:n
                    bin_c = rho_vec(kk);
                    if bin_c == -1 % a => b
                        s_tmp = 'm';
                    elseif bin_c == 0% a => a or b => b 
                        if c_vec(kk) == 1 % a
                            s_tmp = 'a';
                        elseif c_vec(kk) == 0 % b
                            s_tmp = 'b';
                        end
                    elseif bin_c == 1 % b => a
                        s_tmp = 'p';
                    end
                    rho_tmp = strcat(rho_tmp,s_tmp);
                end
            rho(ii,jj) = sym(rho_tmp); 
            end
        end
    end
    % rho = rho_box(n)
      
    end % methods (Static)
end % classdef