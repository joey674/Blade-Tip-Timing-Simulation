% MDOF_Bound function: return the weight area of the peaks and the
% boundary area of the peaks. the function will first oversmooth the data
% to minimize the noise; then with a really low prominence to find peaks so
% that no peaks will be missed; after comparing with peaks idx with
% input(peaks will be little moving because of smoothing the data) we can
% find the coresponded peaks and use the peaks info to find boundary.

% Input variables:
%   - magn: 
%   - peaks_idx: 

% output variables:
%   - weights:
%   - boundary:

function boundary_idx = MDOF_Bound(magn,selected_peaks_idx)
    selected_peaks_idx = sort(selected_peaks_idx);
    
    % get peaks info 
    magn_fp = smoothdata(magn,'gaussian',2000);
    magn_fp = magn_fp * max(magn)/max(magn_fp);% 高度平滑也不能影响总体scale
    min_peak_prominence = 0.05 * (max(magn_fp) - mean(magn_fp));
    [~, locs, widths, ~] = findpeaks(magn_fp, 'MinPeakProminence', min_peak_prominence);
    
    % boundary
    [min_distance, closest_idx] = min(abs(locs - selected_peaks_idx(1)));
    if min_distance > 1000 %如果离太远就有问题了，使用最左或最右
        fprintf("error! select peak no found\n");
        left_boundary_idx = 1;
    else
        peak_loc = locs(closest_idx);
        peak_width = widths(closest_idx);   
        left_boundary_idx = max(peak_loc - 2 * peak_width, 1);% 向左搜索峰的边界
        left_boundary_idx = round(left_boundary_idx);
        [~, minIdx] = min(magn_fp(left_boundary_idx:peak_loc));
        left_boundary_idx = left_boundary_idx + minIdx - 1; 
    end
    clear peak_loc peak_width closest_idx min_distance

    [min_distance, closest_idx] = min(abs(locs - selected_peaks_idx(end)));
    if min_distance > 1000 %如果离太远就有问题了，使用最左或最右
        fprintf("error! select peak no found\n");
        right_boundary_idx = length(magn);
    else
        peak_loc = locs(closest_idx(end));
        peak_width = widths(closest_idx(end));
        right_boundary_idx = min(peak_loc + 2 * peak_width, length(magn_fp));% 向右搜索峰的边界
        right_boundary_idx = round(right_boundary_idx);
        [~, minIdx] = min(magn_fp(peak_loc:right_boundary_idx));
        right_boundary_idx = peak_loc + minIdx - 1;    
    end
    boundary_idx = [round(left_boundary_idx),round(right_boundary_idx)];
end