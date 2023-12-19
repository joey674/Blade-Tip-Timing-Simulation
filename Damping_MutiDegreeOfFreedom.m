function Damping_MutiDegreeOfFreedom(dataset_path,blade,n_modes,data_smoothness,peak_interval_width)
fprintf('[**********damping:muti degree of freedom approximation starts.**********]\n');

n_blades = length(blade);

%% get peak intervals
peak_intervals = MDOF_FindPeakIntervals(n_blades,n_modes,blade,data_smoothness,peak_interval_width);
fprintf('%d ',peak_intervals.');
fprintf('\n');


%% start algorithm
% for saving to file
freq_sorted_tofile = cell(1, n_blades);
magn_sorted_tofile = cell(1, n_blades);
err_sorted_tofile = cell(1, n_blades);
freq_used_tofile = cell(1, n_blades);
params_fitted_tofile = cell(1, n_blades);
model_fitted_tofile = cell(1, n_blades);
% for plot phase
phase_arr = cell(n_modes, 1);
% start iterate(find best fitted freq range)
for blade_idx = 1:n_blades
% for blade_idx = 3:2:9
    fprintf('blade:%d\n',blade_idx);
    % import data;
    blade_data = blade{blade_idx};
    freq = [blade_data.freq];
    magn = [blade_data.magn]; 
    err  = [blade_data.err]; 

    % smooth data 
    magn = smoothdata(magn,'movmean',data_smoothness);
    err = smoothdata(err,'movmean',data_smoothness);

    % find peaks
    peaks_y = [];
    peaks_idx = [];    
    % pks保存的是峰值;locs保存的是峰值对应的magn的idx
    peak_prominence = 0.2 * (max(magn) - mean(err));% 多突出的峰值才被认定为峰值;越小找出的peak越多,越可能获取到peak.但是有可能获取到太多从而干扰
    [pks, locs] = findpeaks(magn, 'MinPeakProminence', peak_prominence);   
    for i = 1:size(peak_intervals, 1)
        peak_interval = peak_intervals(i, :);
        locs_idx = find(freq(locs) >= peak_interval(1) & freq(locs) <= peak_interval(2));% loc_idx就表示的是第几个pks/locs
        % check if we have a peak in the interval
        if ~isempty(locs_idx)
            % select the max peak
            [peak_value_in_interval, locs_idx_idx] = max(pks(locs_idx));
            max_idx = locs_idx(locs_idx_idx);  
            % save the peak
            peak_threshold = 1.4 * mean(err);
            if peak_value_in_interval > peak_threshold && peak_value_in_interval > mean(magn)% if peak is not higher than err, then dont use this peak
                peaks_y(end+1) = pks(max_idx);
                peaks_idx(end+1) = locs(max_idx);
            end
        end
    end
    
    % remove unneed part before first peak and after last peak(ignored_ratio)
    ignored_ratio_best = 0;
    residual_sum_best = Inf; 
    %find ignored_ratio range and step
    ignored_ratio_min = 0.3*peak_threshold / max(peaks_y);
    ignored_ratio_max = max(peak_threshold, 0.3*min(peaks_y)) / max(peaks_y);% 取mean(err)或者min(peaks)的较小的值,但保证大于ignored_ratio_min
    ignored_ratio_step = (ignored_ratio_max-ignored_ratio_min)/10;% will iterate 10 time
    fprintf('min:%d,max:%d,step:%d\n',ignored_ratio_min,ignored_ratio_max,ignored_ratio_step);
    for ignored_ratio_tmp = ignored_ratio_min:ignored_ratio_step:ignored_ratio_max
        freq_tmp = freq;
        magn_tmp = magn;
        err_tmp = err;
        peaks_idx_tmp = peaks_idx;
        if ~isempty(peaks_idx_tmp)
            % ignore the part that lower than ignore_percentage of max peak before first peak
            before_first_peak_idxs = find(freq_tmp < freq_tmp(peaks_idx_tmp(1)));
            threshold_idx = find(magn_tmp(before_first_peak_idxs) < ignored_ratio_tmp * max(peaks_y), 1, 'last');       
            if ~isempty(threshold_idx)
                valid_idx = before_first_peak_idxs(threshold_idx):length(freq_tmp);
                freq_tmp = freq_tmp(valid_idx);
                magn_tmp = magn_tmp(valid_idx);
                err_tmp = err_tmp(valid_idx);     
                peaks_idx_tmp = peaks_idx_tmp - valid_idx(1) + 1;
            end
            % ignore the part that lower than ignore_percentage of max peak after last peak
            after_last_peak_idxs = find(freq_tmp > freq_tmp(peaks_idx_tmp(end)));
            threshold_idx = find(magn_tmp(after_last_peak_idxs) < ignored_ratio_tmp * max(peaks_y), 1, 'first');      
            if ~isempty(threshold_idx)
                valid_idx = 1:(after_last_peak_idxs(threshold_idx) - 1);
                freq_tmp = freq_tmp(valid_idx);
                magn_tmp = magn_tmp(valid_idx);
                err_tmp = err_tmp(valid_idx);     
                peaks_idx_tmp = peaks_idx_tmp(peaks_idx_tmp <= valid_idx(end));
            end
        end        
        params_fitted = MDOF_LMAlgorithm(freq_tmp,magn_tmp,peaks_idx_tmp);
        residual_sum = sum(abs(MDOF_ModelMagn(params_fitted,freq_tmp) - magn_tmp))/length(freq_tmp);
        % fprintf('ignored_ratio_residual_sum:%d\n',residual_sum);
        % save result if this ignored_ratio is better
        if residual_sum < residual_sum_best
            residual_sum_best = residual_sum;
            ignored_ratio_best = ignored_ratio_tmp;   
            freq_best = freq_tmp;
            magn_best = magn_tmp;
            params_fitted_best = params_fitted;       
            peaks_idx_best = peaks_idx_tmp;
        end                   
    end
    fprintf('ignored_ratio_best:%d\n',ignored_ratio_best);
   
    % plot
    figure('units','normalized','outerposition',[0 0 0.7 0.7]);  
    set(gcf, 'WindowStyle', 'docked');
    subplot(2, 1, 1);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (mm)');
    hold on; 
    plot(freq, magn, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Magnitude');
    plot(freq, err, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Error');
    plot(freq(peaks_idx), peaks_y, 'bo', 'DisplayName', 'Peaks');
    legend;       
    plot(freq_best, MDOF_ModelMagn(params_fitted_best,freq_best), 'g--', 'DisplayName', 'Fitted Model');
    hold off;
    subplot(2, 1, 2); 
    plot(freq_best, MDOF_ModelMagn(params_fitted,freq_best) - magn_best, 'm-','DisplayName', 'Residual');
    legend;
    xlabel('Frequency (Hz)');
    ylabel('Residual');
    hold off;
    
    % save to file
    freq_sorted_tofile{blade_idx} = freq;
    magn_sorted_tofile{blade_idx} = magn;
    err_sorted_tofile{blade_idx} = err;
    freq_used_tofile{blade_idx} = freq_best;
    params_fitted_tofile{blade_idx} = reshape(params_fitted_best,4, length(peaks_idx_best)).';
    model_fitted_tofile{blade_idx} = MDOF_ModelMagn(params_fitted_best,freq_best);

    % plot phase 
    phase = MDOF_ModelPhase(params_fitted_best,freq_best);
    % figure('units','normalized','outerposition',[0 0 0.7 0.7]);  
    % set(gcf, 'WindowStyle', 'docked');
    % plot(freq_best, phase, 'g--', 'DisplayName', 'Fitted Model');
    
    for idx = 1:length(peaks_idx_best)
        for i = 1:size(peak_intervals, 1)
            peak_interval = peak_intervals(i, :);
            if (freq_best(peaks_idx_best(idx))>= peak_interval(1)) && (freq_best(peaks_idx_best(idx))<= peak_interval(2))
                phase_arr{i} = [phase_arr{i}, phase(peaks_idx_best(idx))];
            end
        end
    end
end
filename = strrep(dataset_path, '.mat', '_MDOF.mat');
save(filename, 'freq_sorted_tofile','magn_sorted_tofile','err_sorted_tofile','freq_used_tofile','params_fitted_tofile','model_fitted_tofile');

% plot phase by modes
% for i = 1:n_modes   
%     figure('units','normalized','outerposition',[0 0 0.7 0.7]); 
%     set(gcf, 'WindowStyle', 'docked');
%     polarplot(phase_arr{i},2, 'o');
%     hold off;
% end

fprintf('[**********damping:muti degree of freedom approximation finished.**********]\n');
end

