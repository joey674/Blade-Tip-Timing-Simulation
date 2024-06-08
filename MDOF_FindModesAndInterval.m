function [n_modes,peak_intervals] = MDOF_FindModesAndInterval(blade)
    n_blades = length(blade);
    for i = 1:n_blades
        lens(i) = length([blade{i}.magn]);
    end
    magn_all = NaN(max(lens),n_blades);
    for i = 1:n_blades
        magn_all(1:length([blade{i}.magn]),i) = [blade{i}.magn]./max([blade{i}.magn]);
    end
    magn_max = max(magn_all,[],2,"omitnan");
    magn_max = smoothdata(magn_max,'gaussian',1000); % 高度平滑 用来查找适中的interval 
    
    % set peak_prominence 
    freq = [blade{1}.freq];
    min_peak_prominence = 0.1 * (max(magn_max) - mean(magn_max));    

    % find peak 
    [pks, locs,widths,prominences] = findpeaks(magn_max, MinPeakProminence = min_peak_prominence);             
    % fprintf('widths'); fprintf(' %d ',widths); fprintf('\n');
    % fprintf('prominences'); fprintf(' %d ',prominences); fprintf('\n');
    valid_idx = pks > 0.2*(max(magn_max)-mean(magn_max)) + mean(magn_max);
    pks = pks(valid_idx);
    locs = locs(valid_idx);

    % calculate n_modes
    n_modes = length(pks);   

    % set peak_interval_width   
    if length(locs) > 1
        peak_distances = diff(freq(locs));
        min_distance = min(peak_distances);
        peak_interval_width = min_distance;
    else
        fprintf("it is SDOF data.");
        peak_interval_width = 50; 
    end
    
    %% calculate peak interval
    peak_intervals = zeros(n_modes,2);
    for i = 1:n_modes
        lower_bound = blade{1}(locs(i)).freq - peak_interval_width/2;
        upper_bound = blade{1}(locs(i)).freq + peak_interval_width/2;
        peak_intervals(i, :) = [lower_bound, upper_bound];
    end
    
    %% plot
    figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
    set(gcf, 'WindowStyle', 'docked');
    hold on;
    % 绘制幅值数据
    plot(freq, magn_max, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Magnitude');
    % 绘制峰值
    plot(freq(locs), pks, 'ro', 'DisplayName', 'Peaks');
    % 绘制峰值区间
    for i = 1:n_modes
        fill([peak_intervals(i, 1), peak_intervals(i, 2), peak_intervals(i, 2), peak_intervals(i, 1)], ...
             [min(magn_max), min(magn_max), max(magn_max), max(magn_max)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Peak Interval');
    end 
    % 显示图例
    legend('show');
    hold off;   
end



