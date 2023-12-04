function peak_intervals = find_peak_intervals(n_blades, n_modes, blade)
    peak_prominence_max = 20;
    peak_prominence_step = 0.5;
    peak_interval_width = 10; 

    peak_freqs = []; 
    for i = 1:n_blades
        % inport a blade
        blade_data = blade{i};
        freq = [blade_data.freq];
        magn = [blade_data.magn];        
        % smooth 
        magn = smoothdata(magn, 'movmean', 100);  
        % init PeakProminence       
        pks = [];
        locs = [];   
        % try to find min PeakProminence 
        peak_prominence = peak_prominence_max;
        while length(pks) < n_modes && peak_prominence > 0
            [pks, locs] = findpeaks(magn, 'MinPeakProminence', peak_prominence);
            peak_prominence = peak_prominence - peak_prominence_step; % 减小 MinPeakProminence
        end   
        % 添加峰值频率到数组
        peak_freqs = [peak_freqs, freq(locs)];
        % fprintf('blade:%d\n',i);
        % fprintf('peak_prominence:%d\n',peak_prominence);
        % fprintf('peak_freqs:%d\n',freq(locs));
    end
    % sort the freq upward
    peak_freqs = sort(peak_freqs);


    % find intervals (moving interval)
    peak_intervals = zeros(2, n_modes);
    for i = 1:n_modes
        % 把peak_freqs一个一个放进前面这个函数里,返回一个数组,数组里是对应的每个函数的输出(这里输出是每个peak_freqs对应的包含peak个数)
        % max函数返回数组最大值以及其索引,这里我们不考虑最值,只要这个索引.这个索引就是peak_freqs的索引
        [~, idx] = max(arrayfun(@(x) sum(peak_freqs >= x & peak_freqs <= x + peak_interval_width), peak_freqs));
        if isempty(idx)
            break;
        end
        lower_bound = peak_freqs(idx);
        upper_bound = lower_bound + peak_interval_width;
        peak_intervals(:, i) = [lower_bound, upper_bound];
        peak_freqs(peak_freqs >= lower_bound & peak_freqs <= upper_bound) = []; % delete freq we already used 
        % fprintf('[%d,%d]\n',lower_bound,upper_bound);
    end
    peak_intervals = sortrows(peak_intervals', 1);% sort the intervals upward
end



