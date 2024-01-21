% MDOF_LMAlgorithm 

% Input variables:
%   - magn: 
%   - peaks_idx: 

% Adjust variables:
%   - weights_value: 
%   - core_weights_value: 

function params_fitted = MDOF_LMAlgorithm(freq,magn,peaks_idx,weights_idx) 
    % calculate significant mode number for this balde
    n_modes = length(peaks_idx);
    % weights
    weights = ones(size(freq));
    weights_value = 20;
    for i = 1:size(weights_idx, 1)
        start_idx = max(1, weights_idx(i, 1));
        end_idx = min(length(freq), weights_idx(i, 2));
        weights(start_idx:end_idx) = weights_value;
    end
    % set initial value for params (W_m: exci_freq; D_m: damping ratio; r_re; r_im;)
    params_initial = zeros(1, n_modes*4);
    lb = [];ub = [];
    for i = 1:n_modes 
        params_initial((i-1)*4 + 1) = freq(peaks_idx(i));
        params_initial((i-1)*4 + 2) = 0.0001; 
        % r_re does matters and must set to magn_each_blade_smoothed(peaks_locs(i)),r_im is not that important
        params_initial((i-1)*4 + 3) = magn(peaks_idx(i));
        params_initial((i-1)*4 + 4) = 100;
        lb_s = [freq(peaks_idx(i))-5, 0, -Inf, -Inf];
        ub_s = [freq(peaks_idx(i))+5, 0.01, Inf, Inf];
        lb = [lb,lb_s];
        ub = [ub,ub_s];
    end    

    % set lsqnonlin use levenberg-marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
    residual = @(P) weights .* (MDOF_Model(P,freq) - magn);
    params_fitted = lsqnonlin(residual, params_initial, lb, ub, options);   
end

