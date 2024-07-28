function params_fitted = LMAlgorithmSDOF(freq,magn,peaks_idx,weights_idx) 
    %% weights
    weights = ones(size(freq));
    peak_weight = 150;
    edge_weight = 10;
    for i = 1:size(weights_idx, 1)
        peak_idx = peaks_idx(i);
        start_idx = max(1, weights_idx(i, 1));
        end_idx = min(length(freq), weights_idx(i, 2));
        % Gaussian weight distribution from peak to edges
        sigma = (end_idx - start_idx) / 10; 
        x = start_idx:end_idx;
        gaussian_weights = peak_weight * exp(-((x - peak_idx).^2 / (2 * sigma^2)));
        % Ensure edge weights are not less than edge_weight
        gaussian_weights(gaussian_weights < edge_weight) = edge_weight;
        weights(start_idx:end_idx) = gaussian_weights;
    end

    %% set initial value for params (W_m: exci_freq; D_m: damping ratio; K_j:weight factor )
    deviation_peaks_idx = 10;
    params_initial = zeros(1, 3);
    lb = [];ub = [];

    params_initial(1) = freq(peaks_idx);
    params_initial(2) = 0.0001; 
    params_initial(3) = max(magn) * 100000;
    lb_s = [freq(peaks_idx)-deviation_peaks_idx, 0, -Inf];
    ub_s = [freq(peaks_idx)+deviation_peaks_idx, 0.01, Inf];
    lb = [lb,lb_s];
    ub = [ub,ub_s];  

    % set lsqnonlin use levenberg-marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
        'MaxIterations', 1000,'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6);
    residual = @(P) weights .* ( abs(ModelSDOF(P,freq)) - magn ).^2;
    params_fitted = lsqnonlin(residual, params_initial, lb, ub, options);   
end

