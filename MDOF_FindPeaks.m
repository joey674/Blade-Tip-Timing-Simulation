function peaks_idx = MDOF_FindPeaks(magn,freq,err,peak_intervals,n_modes)
    peaks_idx = [];    

    min_peak_prominence = 0.1 * (max(magn) - mean(err));% 初始prominence值(设置稍微偏小)
    [pks, locs] = findpeaks(magn, 'MinPeakProminence', min_peak_prominence); 

    % 对剩下的peak进行是否在interval内的筛选,且只在interval内选择最大的那个
    for i = 1:size(peak_intervals, 1)
        peak_interval = peak_intervals(i, :);
        locs_idx = find(freq(locs) >= peak_interval(1) & freq(locs) <= peak_interval(2));% loc_idx就表示的是第几个pks/locs
        % check if we have a peak in the interval
        if ~isempty(locs_idx)
            % select the max peak
            [peak_value_in_interval, locs_idx_idx] = max(pks(locs_idx));
            max_idx = locs_idx(locs_idx_idx);  
            % 这里设置不仅要大于err的条件,还要是突出的部分最起码要在平均值以上
            if peak_value_in_interval > max( 1.3*mean(err), 0.1*(max(magn)-mean(magn))+mean(magn) )
                peaks_idx(end+1) = locs(max_idx);
            end
        end
    end

    % 把peak_idx从小到大排序 并删除重复的值
    peaks_idx = unique(peaks_idx);
    [peaks_idx, ~] = sort(peaks_idx);
end