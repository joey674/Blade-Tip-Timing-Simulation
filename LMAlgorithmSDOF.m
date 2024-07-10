function params_fitted = LMAlgorithmSDOF(freq,magn,peaks_idx,weights_idx) 
    %% set initial value for params (W_m: exci_freq; D_m: damping ratio; K_j:weight factor )
    deviation_peaks_idx = 10;
    params_initial = zeros(1, 3);
    lb = [];ub = [];

    params_initial(1) = freq(peaks_idx);
    params_initial(2) = 0.0001; 
    params_initial(3) = max(magn) * 400000;
    lb_s = [freq(peaks_idx)-deviation_peaks_idx, 0, -Inf];
    ub_s = [freq(peaks_idx)+deviation_peaks_idx, 0.01, Inf];
    lb = [lb,lb_s];
    ub = [ub,ub_s];  

    % set lsqnonlin use levenberg-marquardt algorithm
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', ...
        'MaxIterations', 1000,'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6);
    residual = @(P) ( abs(ModelSDOF(P,freq)) - magn ).^2;
    params_fitted = lsqnonlin(residual, params_initial, lb, ub, options);   
end

