        % obj = xyz2pmz(obj);
        function obj = xyz2pmz(obj);
            % obj = xyz2pmz(obj)
            % conversion from Cartesian operator basis to lowring/raising operator basis
            axis_in = obj.axis;
            coef_in = obj.coef;
            Ncoef_in = obj.Ncoef;
            spin_no = size(axis_in,2);
        
            axis_out = [];
            coef_out = [];
            for ii = 1:size(axis_in,1)
                axis_tmp = axis_in(ii,:);
        
                % Conversion from Ix and Iy to Ip and Im
                xn = length(find(axis_tmp == 1)); % Example: IxSyKx =>(Ip + Im)(Sp - Sm)(Kp + Km)
                yn = length(find(axis_tmp == 2)); % => 2^2 * 2^1 = 8.
                xyn = 2^(xn + yn);
                axis_out_tmp = repmat(axis_tmp,xyn,1);
                axis_out_tmp(axis_out_tmp == 1) = 0; % Remove Ix
                axis_out_tmp(axis_out_tmp == 2) = 0; % Remove Iy
                coef_out_tmp = ones(xyn,1); %
        
                dec = [0:xyn-1];
                bin = dec2bin(dec,(xn + yn));%str            
                % Converting from str to double for each character
                bin_mat = zeros(xyn,xn + yn);
                for nn = 1:xyn
                    for kk = 1:(xn + yn)
                        bin_mat(nn,kk) = str2double(bin(nn,kk));
                    end
                end
                bin_mat (bin_mat == 0) = 4;
                bin_mat (bin_mat == 1) = 5;

                int_count = 0;
                for jj = 1:spin_no
                    axis_v = axis_tmp(jj);
                    if axis_v == 1 || axis_v == 2
                        int_count = int_count + 1;
                        bin_vec = bin_mat(:,int_count);
                        axis_out_tmp(:,jj) = bin_vec;

                        c_tmp = zeros(size(bin_vec));
                        if axis_v == 1
                            c_tmp(bin_vec == 4) = 1/2;
                            c_tmp(bin_vec == 5) = 1/2;
                        elseif axis_v == 2
                            c_tmp(bin_vec == 4) = 1/(2*1i);
                            c_tmp(bin_vec == 5) = -1/(2*1i);
                        end
                        coef_out_tmp = coef_out_tmp.*c_tmp;
                    end
                end
                axis_out = [axis_out;axis_out_tmp];
                coef_out_tmp = coef_out_tmp*coef_in(ii)*Ncoef_in(ii);
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
        
            obj2 = PO();
            obj2.axis = axis_out;
            obj2.coef = coef_out;
            obj2.spin_label = obj.spin_label;
            obj2.bracket = bracket_out;% change the property of bracket later.
            obj = CombPO(obj2);
            % obj.coef = rewrite(obj.coef,'exp');
        end
        % xyz2pmz
