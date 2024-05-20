function indices_downsampled = custom_downsample(indices, factor)
    % 自定义的下采样函数，接受每100个元素中需要取的元素数量
    num_elements = length(indices); % 元素总数
    group_size = 100; % 每组元素的大小
    num_to_pick = round(factor); % 每组中需要选取的元素数量
    
    % 初始化输出数组
    indices_downsampled = [];
    
    % 遍历所有组
    for start_idx = 1:group_size:num_elements
        end_idx = min(start_idx + group_size - 1, num_elements); % 确保不超出数组界限
        current_group = indices(start_idx:end_idx); % 当前组
        if num_to_pick > length(current_group)
            num_to_pick = length(current_group); % 如果组中的元素不足以取出指定数量，调整数量
        end
        
        % 从当前组中均匀地选取元素
        pick_indices = round(linspace(1, length(current_group), num_to_pick));
        indices_downsampled = [indices_downsampled, current_group(pick_indices)];
    end
end


