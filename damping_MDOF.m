dataset = 'App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat';
n_blades = 12;
n_modes = 4;
EO = 24;
fprintf('[**********damping ratio calculation start.**********]\n');


%% load RPM/Magn/Err from file
data = load(dataset, 'mean_RPM' , 'P_Magn', 'Fit_Error');
if isfield(data,'mean_RPM')
    Freq = data.mean_RPM * EO / 60;    
else
    fprintf('mean_RPM does not exist in the file\n');
end
if isfield(data, 'P_Magn')     
    Magn = cell(1, n_blades);     
    for i = 1:n_blades
        Magn{i} = data.P_Magn{i,24,:};  
    end    
else
    fprintf('P_Magn does not exist in the file\n');
end
if isfield(data,'Fit_Error')
    Err = cell(1, n_blades);  
    for i = 1:n_blades
        Err{i} = data.Fit_Error{i,24,:};  
    end    
else
    fprintf('Fit_Error does not exist in the file\n');
end


%% pre process of data
blade = cell(1, n_blades);   
% find max length
max_length = max([length(Freq), cellfun(@length, Magn)]);   
for i = 1:n_blades
    % assign Magn{i} and Freq to the structure
    blade{i} = repmat(struct('freq', NaN, 'magn', NaN, 'err', NaN ), 1, max_length);       
    len_freq = length(Freq);
    len_magn = length(Magn{i});
    for j = 1:max_length
        if j <= len_freq
            blade{i}(j).freq = Freq(j);
        end
        if j <= len_magn
            blade{i}(j).magn = Magn{i}(j);
            blade{i}(j).err = Err{i}(j);
        end
    end    
    % sort the structure array by freq
    [~, idx_sort] = sort([blade{i}.freq]);
    blade{i} = blade{i}(idx_sort);   
    % remove elements contain nan
    idx_remove = isnan([blade{i}.freq]) | isnan([blade{i}.magn]);                   
    blade{i}(idx_remove) = [];   
end


%% get peak intervals
% peak_intervals = [13739, 13748; 13773, 13779; 13800, 13807; 13831, 13839];
% peak_intervals = [1.374005e+04, 1.375005e+04; 1.377473e+04, 1.378473e+04; 1.380403e+04, 1.381403e+04; 1.383250e+04, 1.384250e+04];
peak_intervals = find_peak_intervals(n_blades,n_modes,blade);
fprintf('%d ',peak_intervals);
fprintf('\n');


%% start algorithm
ignored_ratio_max = 0.3;
ignored_ratio_step = 0.05;
peak_prominence = 5;% 多突出的峰值才被认定为峰值;越小找出的peak越多,越可能获取到peak.但是有可能获取到太多从而干扰
peak_smoothness = 100; % 对数据平滑的程度 用来寻找峰值
for blade_idx = 1:n_blades
% for blade_idx = 7:7
    fprintf('blade:%d\n',blade_idx);
    % import data;
    blade_data = blade{blade_idx};
    freq = [blade_data.freq];
    magn = [blade_data.magn]; 
    err  = [blade_data.err]; 

    % smooth data 
    magn = smoothdata(magn,'movmean',peak_smoothness);
    err = smoothdata(err,'movmean',peak_smoothness);

    % find peak
    peaks_y = [];
    peaks_idx = [];    
    % pks保存的是峰值;locs保存的是峰值对应的magn的idx
    [pks, locs] = findpeaks(magn, 'MinPeakProminence', peak_prominence);   
    for i = 1:size(peak_intervals, 1)
        peak_interval = peak_intervals(i, :);
        locs_idx = find(freq(locs) >= peak_interval(1) & freq(locs) <= peak_interval(2));% loc_idx就表示的是第几个pks/locs
        % check if we have a peak in the interval
        if ~isempty(locs_idx)
            % select the max peak
            [~, locs_idx_idx] = max(pks(locs_idx));
            max_idx = locs_idx(locs_idx_idx);  
            % save the peak
            peaks_y(end+1) = pks(max_idx);
            peaks_idx(end+1) = locs(max_idx);
        end
    end

      
    % remove unneed part before first peak and after last peak
    ignored_ratio_best = 0;
    residual_sum_best = Inf;   
    for ignored_ratio_tmp = 0:ignored_ratio_step:ignored_ratio_max
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
        
        params_fitted = LM_Algorithm(freq_tmp,magn_tmp,peaks_idx_tmp);
        residual_sum = sum(abs(LM_Residual(params_fitted,freq_tmp,magn_tmp)))/length(freq_tmp);
        fprintf('ignored_ratio_residual_sum:%d\n',residual_sum);

        % save result if this ignored_ratio is better
        if residual_sum < residual_sum_best
            residual_sum_best = residual_sum;
            ignored_ratio_best = ignored_ratio_tmp;   
            freq_best = freq_tmp;
            magn_best = magn_tmp;
            params_fitted_best = params_fitted;           
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
    plot(freq_best, MDOF_Model(params_fitted_best,freq_best), 'g--', 'DisplayName', 'Fitted Model');
    hold off;
    subplot(2, 1, 2); 
    plot(freq_best, LM_Residual(params_fitted_best,freq_best,magn_best), 'm-','DisplayName', 'Residual');
    legend;
    xlabel('Frequency (Hz)');
    ylabel('Residual');
    hold off;
end

fprintf('[**********damping ratio calculation finished.**********]\n');



