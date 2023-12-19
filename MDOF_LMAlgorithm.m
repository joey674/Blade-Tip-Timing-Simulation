function params_fitted = MDOF_LMAlgorithm(freq,magn,peaks_idx) 
    % calculate significant mode number for this balde
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
    lb_group = [0, 0, -Inf, -Inf];
    ub_group = [Inf, 1, Inf, Inf];
    lb = repmat(lb_group, 1, n_modes);
    ub = repmat(ub_group, 1, n_modes);
    % create weight func
    weights = ones(size(freq)); % 初始化权重为 1
    peak_weight = 10; % 峰值附近的权重

    for i = 1:length(peaks_idx)
        peak_loc = peaks_idx(i);
        % 增加峰值附近的权重
        weights(max(1, peak_loc-50):min(end, peak_loc+50)) = peak_weight;
    end
    

    % set lsqnonlin use levenberg-marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
    residual = @(P) weights .* (MDOF_ModelMagn(P,freq) - magn);
    params_fitted = lsqnonlin(residual, params_initial, lb, ub, options);   
end

