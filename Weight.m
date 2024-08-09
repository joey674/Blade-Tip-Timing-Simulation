function weights_idx = Weight(magn, selected_peaks_idx)
    % find minimal
    magn_oversmoothed = smoothdata(magn, 'gaussian', 2000);
    [~, loc_idx] = findpeaks(-magn_oversmoothed);
    weights_idx = zeros(length(selected_peaks_idx), 2);
    
    % find the nearlest minimals for both sides of every peaks
    for i = 1:length(selected_peaks_idx)
        peak_idx = selected_peaks_idx(i);

        left_loc_idx = loc_idx(loc_idx < peak_idx);
        if ~isempty(left_loc_idx)
            left_valley = left_loc_idx(end); 
        else
            left_valley = 1; 
        end

        right_valleys = loc_idx(loc_idx > peak_idx);
        if ~isempty(right_valleys)
            right_valley = right_valleys(1); 
        else
            right_valley = length(magn); 
        end

        % check if weight overlap
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
        
        % result
        weights_idx(i, :) = [left_valley, right_valley];
    end
    % disp(weights_idx);
end

