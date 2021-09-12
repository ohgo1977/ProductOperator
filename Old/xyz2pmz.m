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
        xn = 2^length(find(axis_tmp == 1)); % Example: IxSyKx =>(Ip + Im)(Sp - Sm)(Kp + Km)
        yn = 2^length(find(axis_tmp == 2)); % => 2^2 * 2^1 = 8.
        xyn = xn*yn;
        axis_out_tmp = zeros(xyn,spin_no); % 0: additive identity
        coef_out_tmp = ones(xyn,1); % 1: multiplicative identity

        for jj = 1:spin_no
            axis_v = axis_tmp(jj);
            if axis_v == 1 % Ix = 1/2*(Ip + Im), 4: Ip, 5: Im
                v_tmp = repmat([4;5],xyn/2,1);  % [4;5;4;5;4;5;...]
                c_tmp = repmat([1;1]/2,xyn/2,1);% [1;1;1;1;1;1;...]/2

            elseif axis_v == 2 % Iy = 1/(2*1i)*(Ip - Im)
                v_tmp = [repmat(4,xyn/2,1);repmat(5,xyn/2,1)];% [4;4;4;...;5;5;5;...]
                c_tmp = [repmat(1/(2*1i),xyn/2,1);repmat(-1/(2*1i),xyn/2,1)];%[1;1;1;...;-1;-1;-1;...]/(2*1i)

            elseif axis_v == 3 % Iz
                v_tmp = repmat(3,xyn,1);
                c_tmp = ones(xyn,1);

            elseif axis_v == 0 
                v_tmp = zeros(xyn,1);
                c_tmp = ones(xyn,1);
            end
            axis_out_tmp(:,jj) = v_tmp;
            coef_out_tmp = coef_out_tmp.*c_tmp;
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
    obj2.coherence = PO.rho_box(spin_no);
    obj2.bracket = bracket_out;% change the property of bracket later.
    obj = CombPO(obj2);
    % obj.coef = rewrite(obj.coef,'exp');
end