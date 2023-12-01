function params_fitted = LM_Algorithm(freq,magn,peaks_idx) 
    % calculate mode number
   n_modes = length(peaks_idx);
    % set initial value for params (W_m: exci_freq; D_m: damping ratio; r_re; r_im;)
    params_initial = zeros(1, n_modes*4);
    for i = 1:n_modes 
        params_initial((i-1)*4 + 1) = freq(peaks_idx(i));
        params_initial((i-1)*4 + 2) = 0.0001; 
        % r_re does matters and must set to magn_each_blade_smoothed(peaks_locs(i)),r_im is not that important
        params_initial((i-1)*4 + 3) = magn(peaks_idx(i));
        params_initial((i-1)*4 + 4) = 1000;
    end
    % set boundary
    lb = [
        0,0,-Inf,-Inf,...
        0,0,-Inf,-Inf,...
        0,0,-Inf,-Inf,...
        0,0,-Inf,-Inf,
        ];
    ub = [
        Inf,1,Inf,Inf,...
        Inf,1,Inf,Inf,...
        Inf,1,Inf,Inf,...
        Inf,1,Inf,Inf,                                        
        ];
    % set lsqnonlin use levenberg-marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');
    residual = @(P) MDOF_Model(P,freq) - magn;
    params_fitted = lsqnonlin(residual, params_initial, lb, ub, options);   
end

