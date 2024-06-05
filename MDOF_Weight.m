% MDOF_Weight function: return the weight area of the peaks and the
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

function weights_idx = MDOF_Weight(magn,selected_peaks_idx)

    %%      
    


    %%
    selected_peaks_idx = sort(selected_peaks_idx);
    weights_idx = [];
    
    % get peaks info 
    magn_fp = magn;
    min_peak_prominence = 0.1 * (max(magn_fp) - mean(magn_fp));%% 找到尽可能多的peak
    [~, locs, widths, ~] = findpeaks(magn_fp, 'MinPeakProminence', min_peak_prominence);
    
    % get through info
    magn_inv = -magn;
    min_peak_prominence_inv = 0.5 * (max(magn_fp) - mean(magn_fp));
    [~, troughs_locs, ~, troughs_prominences] = findpeaks(magn_inv,'MinPeakProminence', min_peak_prominence_inv);
    
    % peaks weight 
    for i = 1:length(selected_peaks_idx)
        [min_distance, closest_idx] = min(abs(locs - selected_peaks_idx(i)));% 找到最近的peaks
        if min_distance > 1000 %如果离太远就有问题了，算法跳过这个峰值并对这个峰值使用默认值
            fprintf("weights: select peak no found\n");
            left_boundary_idx = selected_peaks_idx(i)-100;
            right_boundary_idx = selected_peaks_idx(i)+100;
        else
            peak_loc = locs(closest_idx);
            peak_width = widths(closest_idx);
              
            left_search_range_idx = round(max(peak_loc - 2 * peak_width, 1));% 搜索范围
            right_search_range_idx = round(min(peak_loc + 2 * peak_width, length(magn)));
            left_troughs = find(troughs_locs >= left_search_range_idx & troughs_locs <= peak_loc); % 在搜索范围内找到波谷的索引
            right_troughs = find(troughs_locs >= peak_loc & troughs_locs <= right_search_range_idx);    
            if ~isempty(left_troughs) % 找到左侧和右侧最明显的波谷
                [~, max_prom_idx] = max(troughs_prominences(left_troughs));
                left_boundary_idx = troughs_locs(left_troughs(max_prom_idx));
            else
                left_boundary_idx = left_search_range_idx;
            end
            if ~isempty(right_troughs) 
                [~, max_prom_idx] = max(troughs_prominences(right_troughs));
                right_boundary_idx = troughs_locs(right_troughs(max_prom_idx));
            else
                right_boundary_idx = right_search_range_idx;
            end
        end

        % 对获取的左右边界进行检测:是否与别的组重合或者越过了别的peaks
        if i ~= length(selected_peaks_idx)% 右边peaks[0,3/4]
            if right_boundary_idx > selected_peaks_idx(i)+(selected_peaks_idx(i+1) - selected_peaks_idx(i))*3/4
                right_boundary_idx = selected_peaks_idx(i)+(selected_peaks_idx(i+1) - selected_peaks_idx(i))*1/2;
                fprintf("weights: right_boundary_idx exceed\n");
            end
        end
        if i ~= 1
            if left_boundary_idx < weights_idx(end,2)
                left_boundary_idx = weights_idx(end,2);
                fprintf("weights: left_boundary_idx exceed\n");
            end
        end
        weights_idx = [weights_idx;[round(left_boundary_idx),round(right_boundary_idx)]];   
    end
    clear peak_loc peak_width left_boundary_idx right_boundary_idx

    %  % peaks weight 
    % for i = 1:length(selected_peaks_idx)
    %     [min_distance, closest_idx] = min(abs(locs - selected_peaks_idx(i)));% 找到最近的peaks
    %     if min_distance > 1000 %如果离太远就有问题了，算法跳过这个峰值并对这个峰值使用默认值
    %         fprintf("error! select peak no found\n");
    %         weights_idx = [weights_idx;[selected_peaks_idx(i)-50,selected_peaks_idx(i)+50]];
    %         left_boundary_idx = selected_peaks_idx(i)-100;
    %         right_boundary_idx = selected_peaks_idx(i)+100;
    %     else
    %         peak_loc = locs(closest_idx);
    %         peak_width = widths(closest_idx);
    % 
    %         left_search_range_idx = max(peak_loc - 2 * peak_width, 1);% 向左搜索峰的边界
    %         left_search_range_idx = round(left_search_range_idx);
    %         [~, min_idx] = min(magn_fp(left_search_range_idx:peak_loc));
    %         left_boundary_idx = max(left_search_range_idx + min_idx - 1, peak_loc - 1/2 * peak_width); %这个boundary是最靠近的凹陷处,或者是峰宽的一半,取最小;            
    % 
    %         right_search_range_idx = min(peak_loc + 2 * peak_width, length(magn_fp));% 向右搜索峰的边界
    %         right_search_range_idx = round(right_search_range_idx);
    %         [~, min_idx] = min(magn_fp(peak_loc:right_search_range_idx));
    %         right_boundary_idx = min(peak_loc + min_idx - 1,peak_loc + 1/2 * peak_width); %这个boundary是最靠近的凹陷处,或者是峰宽的一半,取最小;            
    %     end
    %     % 对获取的左右边界进行检测:是否与别的组重合或者越过了别的peaks
    %     if i ~= length(selected_peaks_idx)
    %         if right_boundary_idx > selected_peaks_idx(i+1)
    %             right_boundary_idx = (selected_peaks_idx(i) + selected_peaks_idx(i+1))/2;
    %             fprintf("error! right_boundary_idx exceed\n");
    %         end
    %     end
    %     if i ~= 1
    %         if left_boundary_idx < weights_idx(end,2)
    %             left_boundary_idx = weights_idx(end,2);
    %             fprintf("error! left_boundary_idx exceed\n");
    %         end
    %     end
    % 
    %     % 
    %     weights_idx = [weights_idx;[round(left_boundary_idx),round(right_boundary_idx)]];   
    % end
    % clear peak_loc peak_width left_boundary_idx right_boundary_idx
end