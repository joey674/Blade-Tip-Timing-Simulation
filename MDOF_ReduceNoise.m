function magn_smoothed = MDOF_ReduceNoise(magn, freq, err)
    %% 查找密度范围
    density_levels = [1, 0.4, 0.1, 0.05, 0.01]; 
    downsample_factors = [10, 15, 20, 40, 50]; 
    [density_err, values_err] = ksdensity(err);

    [density_max, index_max] = max(density_err);
    value_lower_last_err = values_err(index_max);
    value_upper_last_err = values_err(index_max);
    freq_downsampled = freq;
    magn_downsampled = magn;      
    for i = 1:length(density_levels)
        target_density = density_levels(i) * density_max;
        index_left = index_max;
        while index_left > 1 && density_err(index_left) >= target_density
            index_left = index_left - 1;
        end
        index_right = index_max;
        while index_right < length(density_err) && density_err(index_right) >= target_density
            index_right = index_right + 1;
        end
        value_lower_current_err = values_err(index_left);
        value_upper_current_err = values_err(index_right);

        %
        indices = [find(magn_downsampled > value_lower_current_err & magn_downsampled < value_lower_last_err),...
                   find(magn_downsampled < value_upper_current_err & magn_downsampled > value_upper_last_err)];

        %
        % 自定义的下采样函数，接受每100个元素中需要取的元素数量
        elements_num = length(indices); % 元素总数
        size_group = 100; % 每组元素的大小
        to_pick_num = round(downsample_factors(i)); % 每组中需要选取的元素数量
                
        % 初始化输出数组
        indices_downsampled = [];
                
        % 遍历所有组
        for start_idx = 1:size_group:elements_num
            end_idx = min(start_idx + size_group - 1, elements_num); % 确保不超出数组界限
            group_current = indices(start_idx:end_idx); 
            if to_pick_num > length(group_current)
                to_pick_num = length(group_current); % 如果组中的元素不足以取出指定数量，调整数量
            end
                    
            % 从当前组中均匀地选取元素
            pick_indices = round(linspace(1, length(group_current), to_pick_num));
            indices_downsampled = [indices_downsampled, group_current(pick_indices)];
        end

        %
        indices_outside = [find(magn_downsampled < value_lower_current_err),...
                           find(magn_downsampled > value_upper_current_err),...
                           find(magn_downsampled > value_lower_last_err & magn_downsampled < value_upper_last_err)];

        indices_all = [indices_downsampled, indices_outside];
                
        magn_downsampled = magn_downsampled(sort(indices_all)); 
        freq_downsampled = freq_downsampled(sort(indices_all)); 
                
        value_lower_last_err = value_lower_current_err;
        value_upper_last_err = value_upper_current_err;
    end

    %% 先平滑数据
    magn_downsampled = smoothdata(magn_downsampled, 'gaussian', 200);
        
    %% 插值以恢复数据长度
    magn_interpolated = interp1(freq_downsampled, magn_downsampled, freq, 'linear', 'extrap');
    
    magn_smoothed = magn_interpolated;  
end


