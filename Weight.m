function weights_idx = Weight(magn, selected_peaks_idx)
    % 对数据高度平滑以找到波谷
    magn_oversmoothed = smoothdata(magn, 'gaussian', 2000);
    
    % 对-magn_oversmoothed调用findpeaks即可找到波谷
    [~, loc_idx] = findpeaks(-magn_oversmoothed);

    % 初始化weights_idx数组
    weights_idx = zeros(length(selected_peaks_idx), 2);

    % 对每个峰值找到左右最近的波谷
    for i = 1:length(selected_peaks_idx)
        peak_idx = selected_peaks_idx(i);

        % 找到左边最近的波谷
        left_loc_idx = loc_idx(loc_idx < peak_idx);
        if ~isempty(left_loc_idx)
            left_valley = left_loc_idx(end); % 左边最近的波谷
        else
            left_valley = 1; % 如果没有左边波谷，设为1
        end

        % 找到右边最近的波谷
        right_valleys = loc_idx(loc_idx > peak_idx);
        if ~isempty(right_valleys)
            right_valley = right_valleys(1); % 右边最近的波谷
        else
            right_valley = length(magn); % 如果没有右边波谷，设为数据的最后一个点
        end

        % 排除跨界的问题
        if i>1
            if left_valley < selected_peaks_idx(i-1)
                left_valley = selected_peaks_idx(i-1);
            end
        end
        if i<length(selected_peaks_idx)
            if right_valley > selected_peaks_idx(i+1)
                right_valley = selected_peaks_idx(i+1);
            end
        end

        % 记录左右波谷的索引
        weights_idx(i, :) = [left_valley, right_valley];
    end
    disp(weights_idx);
end




