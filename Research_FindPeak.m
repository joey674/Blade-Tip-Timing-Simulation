function Research_FindPeak(blade, EO)
    %{
        init params
    %}
    n_blades = length(blade);

    %{
        find peaks interval
        先计算所有叶片数据合并在一起,然后归一化
        再计算每一个freq点的最大值
        观测整体的模态 这个图会显示出所有可能的模态
        高度平滑 用来查找适中的interval 
    %}
    for i = 1:n_blades
        lens(i) = length([blade{i}.magn]);
    end
    magn_all = NaN(max(lens), n_blades);
    for i = 1:n_blades
        magn_all(1:length([blade{i}.magn]), i) = [blade{i}.magn] ./ max([blade{i}.magn]);
    end
    magn_max = max(magn_all, [], 2, "omitnan");
    magn_max = smoothdata(magn_max, 'gaussian', 2000); 

    %{
        set parameter
        先用blade1的freq作为总共的freq 因为所有freq都一样
        峰值之间的最小距离 = min_peak_interval_width/2 = 4/2,
        峰值的最大浮动宽度 = max_peak_interval_width = 8,
    %}
    freq = [blade{1}.freq]; 
    min_peak_prominence = 0.1 * (max(magn_max) - mean(magn_max)); 
    min_peak_height = 0.1 * (max(magn_max) - mean(magn_max)) + mean(magn_max); 
    min_peak_interval_width = 4; 
    max_peak_interval_width = 8;
    n_peaks = 10;   
    min_peak_prominence_tro = 0.01 * (max(magn_max) - mean(magn_max)); 

    %{
        find peaks 
        这里loc返回的就不是index了 是freq
    %}
    [pks, locs, widths, prominences] = findpeaks(magn_max, freq, ...
        'MinPeakProminence', min_peak_prominence, ...
        'MinPeakHeight', min_peak_height, ...
        'MinPeakDistance', min_peak_interval_width / 2, ...
        'NPeaks', n_peaks...
    );

    %{
        find troughs
        这里loc_tro返回的就不是index了 是freq
    %}
    [tro_pks, tro_locs, tro_widths, tro_prominences] = findpeaks(-magn_max, freq, ...
        'MinPeakProminence', min_peak_prominence_tro ...
    );

    %{
        find peak_interval_width 
        最宽为max_peak_interval_width 最小为两峰值最小间距min_peak_distance
        第一层判断:有没有找到峰值 如果是SDOF就直接设置为max_peak_interval_width
        第二层判断峰值的最小距离 如果太大就设置为max_peak_interval_width
    %}
    if length(locs) > 1 
        peak_distances = diff(locs);
        min_distance = min(peak_distances);
        if min_distance < max_peak_interval_width
            peak_interval_width = min_distance;
        else 
            peak_interval_width = max_peak_interval_width; 
        end
    else
        peak_interval_width = max_peak_interval_width; 
    end

    %{
        calculate n_modes
        给后面使用的参数
    %}
    n_modes = length(pks);

    %{
        calculate peak interval
        找到最优两个最近的波谷作为边界
    %}
    peak_intervals = zeros(n_modes, 2);
    for i = 1:n_modes
        % 找到每个波峰左右最近的波谷
        left_trough = max(tro_locs(tro_locs < locs(i)));
        right_trough = min(tro_locs(tro_locs > locs(i)));
        if isempty(left_trough)
            left_trough = locs(i) - min_peak_interval_width / 2;
        end
        if isempty(right_trough)
            right_trough = locs(i) + min_peak_interval_width / 2;
        end
        peak_intervals(i, :) = [left_trough, right_trough];
    end

    %{
        plot
    %}
    figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
    title(sprintf('General Magnitude plot EO%d', EO));
    xlabel('Frequency ');
    ylabel('Magnitude ');
    set(gcf, 'WindowStyle', 'docked');
    hold on;
    plot(freq, magn_max, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Magnitude');
    plot(locs, pks, 'ro', 'DisplayName', 'Peaks');
    plot(tro_locs, -tro_pks, 'bo', 'DisplayName', 'Troughs');
    for i = 1:n_modes
        left_bound = peak_intervals(i, 1);
        right_bound = peak_intervals(i, 2);
        fill([left_bound, right_bound, right_bound, left_bound], ...
             [min(magn_max), min(magn_max), max(magn_max), max(magn_max)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end 
    for i = 1:n_modes
        text(locs(i), pks(i), ...
             sprintf('H: %.2f\nP: %.2f\nW: %.2f', pks(i), prominences(i), peak_intervals(i, 2) - peak_intervals(i, 1)), ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    legend('show');
    hold off;

    peak_freqs = locs;
    
    %{
        n_modes         所有叶片 总共有多少模态
        peak_intervals  所有模态的峰值处于的可能区间
        peak_freqs      所有模态的峰值 
    %}
    clearvars -except blade EO n_blades n_modes peak_intervals peak_freqs


    %% 对每个叶片执行寻找峰值策略
    peaks_all(n_blades) = struct('peaks', struct('magn', {}, 'freq', {}, 'prominences', {}, 'widths', {}));
    for blade_idx = 1:n_blades
        %{
            init and smooth
            平滑方式先用rlowess 之后修改成写好的自制方法
        %}
        disp(['blade:', num2str(blade_idx)]); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 
        magn = smoothdata(magn, 'rlowess', 200); 
    
        %% 根据所有叶片的相似模态形状寻找峰值 
        %{
            set parameters
            这里由于没有过度平滑 所以这里的findpeak获得的参数prominence和width都没有意义 会被噪音干扰
            这里的参数是自定义的findpeak方法,使用波峰波谷相对高差来算prominence
        %}  
        min_peak_height = 0.2 * (max(magn) - mean(err)) + mean(err);
        min_peak_prominence = 0.01 * (max(magn) - mean(err));
        disp(['Min Peak Height: ', num2str(min_peak_height)]);
        peaks = struct('magn', {}, 'freq', {}, 'prominences', {}, 'widths', {});
    
        %{
            find peak
            % 获取当前模态区间
            % 找到区间内的峰值,先使用matlab的findpeak找通用的峰值而不是找最
            大值 相当于第一层小小的筛选
            % 选择最大的峰值，并找到其左右波谷
            % 判断是否满足突出度 (prominence) 条件:左右都要和最小的波谷进行
            比较,左右两边都满足
            % 如果满足条件，将其作为峰值，否则尝试下一个峰值，最多尝试三次;
            避免左右有大波峰 大波峰上的小波峰的干扰;由于设置是左右两边都满足,
            则在一般情况下能完成任务
            % 更新 peaks 结构体数组，包含 magn, freq, prominences, 和 widths
        %}
        for i = 1:n_modes 
            interval = peak_intervals(i, :);
            disp(['Interval: [', num2str(interval(1)), ', ', num2str(interval(2)), ']']);
            
            in_interval = (freq >= interval(1)) & (freq <= interval(2));
            interval_freq = freq(in_interval);
            interval_magn = magn(in_interval);
            
            [pks, locs] = findpeaks(interval_magn, interval_freq);          
            if isempty(pks)
                disp('No peaks found that meet the specified criteria.');
                continue;
            end
        
            found_peak = false;
            for attempt = 1:3
                [~, max_idx] = max(pks);
                new_peak.magn = pks(max_idx);
                new_peak.freq = locs(max_idx);
                disp(['Max peak at frequency: ', num2str(new_peak.freq), ' with magnitude: ', num2str(new_peak.magn)]);
        
                left_trough_magn = min(interval_magn(interval_freq < new_peak.freq));
                right_trough_magn = min(interval_magn(interval_freq > new_peak.freq));
        
                if (new_peak.magn - left_trough_magn > min_peak_prominence) && ...
                   (new_peak.magn - right_trough_magn > min_peak_prominence) && ...
                   (new_peak.magn > min_peak_height)
                    new_peak.prominences = max(new_peak.magn - left_trough_magn, new_peak.magn - right_trough_magn);
                    left_trough_idx = find(interval_magn == left_trough_magn & interval_freq < new_peak.freq, 1, 'last');
                    right_trough_idx = find(interval_magn == right_trough_magn & interval_freq > new_peak.freq, 1, 'first');
                    new_peak.widths = interval_freq(right_trough_idx) - interval_freq(left_trough_idx);
                    peaks = [peaks, new_peak];
                    found_peak = true;
                    break;
                end
       
                pks(max_idx) = [];
                locs(max_idx) = [];
                if isempty(pks)
                    disp('No peaks found that meet the specified criteria.');
                    break;
                end
            end
            if ~found_peak
                disp('No valid peaks found in this interval.');
            end
        end


    
        %{
            plot
        %}
        figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
        title(sprintf('EO%d blade%d', EO, blade_idx));
        xlabel('Frequency ');
        ylabel('Magnitude ');
        set(gcf, 'WindowStyle', 'docked');
        hold on;
        plot(freq, magn, 'Color', [0.8, 0.9, 1.0], 'DisplayName', 'Magnitude');
        for j = 1:length(peaks)
            plot(peaks(j).freq, peaks(j).magn, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Selected Peak');
            text(peaks(j).freq, peaks(j).magn, sprintf('M: %.2f\nP: %.2f\nW: %.2f', ...
                peaks(j).magn, peaks(j).prominences, peaks(j).widths), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
        hold off;
        
        %{ 
            获得结果数组
        %}
        peaks_all(blade_idx).peaks = peaks;


        %% 使用相位密度寻找峰值
        % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');
        % title(sprintf('EO%d, blade%d', EO, blade_idx));xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
        % hold on;
        % legend;  
        % 
        % subplot(4,1,1);
        % plot(freq, magn,'o','Color', [0.8, 0.9, 1.0]);  hold on;
        % plot(freq, smoothdata(magn,'rlowess',200),'o');  
        % 
        % subplot(4,1,2);
        % plot(freq, phase,'o','Color', [0.8, 0.9, 1.0]);  
        % 
        % % 定义频率范围和区间数
        % num_bins = 50;  %%%%%%%%%%%%%%%%%%%% 区间精度,分得越多得到的范围越精细
        % freq_edges = linspace(min(freq), max(freq), num_bins + 1);
        % phase_range = [-0.2, 0.2]; %%%%%%%%%%%%%%%%% 相位范围
        % 
        % % 初始化密度数组
        % phase_density = zeros(1, num_bins);
        % freq_centers = (freq_edges(1:end-1) + freq_edges(2:end)) / 2;
        % 
        % % 计算每个频率区间内相位在指定范围内的密度
        % for i = 1:num_bins
        %     % 当前频率区间
        %     idx = freq >= freq_edges(i) & freq < freq_edges(i+1);
        %     current_phase = phase(idx);
        % 
        %     % 计算相位在-0.5到0.5之间的点数比例作为密度
        %     num_points_in_range = sum(current_phase >= phase_range(1) & current_phase <= phase_range(2));
        %     total_points = length(current_phase);
        % 
        %     if total_points > 0
        %         phase_density(i) = num_points_in_range / total_points;
        %     else
        %         phase_density(i) = 0;
        %     end
        % end
        % 
        % % % 绘制相位密度图
        % subplot(4,1,3);
        % plot(freq_centers, phase_density, 'o', 'LineWidth', 2);
        % 
        % % 获取所有密度在0.05以下的频率区间的范围
        % low_density_threshold = 0.1;%%%%%%%%%%%%%%%%%%%%%%相位密度,中间空的地方就是有相变
        % low_density_ranges = [];
        % 
        % for i = 1:num_bins
        %     if phase_density(i) < low_density_threshold
        %         low_density_ranges = [low_density_ranges; freq_edges(i), freq_edges(i+1)];
        %     end
        % end
        %% 使用findpeak函数寻找峰值
        % subplot(4,1,4);
        % 
        % % 降噪
        % magn = Damping_NoiseFilter(magn);
        % plot(freq, magn,'Color', [0.8, 0.9, 1.0]); hold on; 
        % 
        % % 
        % min_prominence = 1/4*mean(magn);%%%%%%%%%%%%%%%%%%%%%峰值的显著性 这个值可以小一点,多收入一些peak
        % min_peakdistance = 1/30*(max(freq)-min(freq));%%%%%%%%%%%%%%%%%%峰值之间的距离 这个值可以设置得大一些,避免peak挤在一起
        % [pks, locs] = findpeaks(magn, freq, 'MinPeakProminence', min_prominence,'MinPeakDistance',min_peakdistance);
        % plot(locs, pks, 'o', 'DisplayName', 'Peaks through findpeak');
        % 
        % %
        % % 标注低密度频率区间(峰值所在区间)
        % for i = 1:size(low_density_ranges, 1)
        %     x_rect = [low_density_ranges(i, 1), low_density_ranges(i, 2), low_density_ranges(i, 2), low_density_ranges(i, 1)];
        %     y_rect = [min(magn), min(magn), max(magn), max(magn)];
        %     fill(x_rect, y_rect, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Low Density Region for peak');
        % end
        %% 将两个办法结合在一起筛选
        % % 筛选在低密度范围内的峰值
        % selected_pks = [];
        % selected_locs = [];
        % 
        % for i = 1:length(locs)
        %     for j = 1:size(low_density_ranges, 1)
        %         if locs(i) >= low_density_ranges(j, 1) && locs(i) < low_density_ranges(j, 2)
        %             selected_pks = [selected_pks, pks(i)];
        %             selected_locs = [selected_locs, locs(i)];
        %             break;
        %         end
        %     end
        % end
        % 
        % % 绘制筛选出的峰值
        % plot(selected_locs, selected_pks, 'bo', 'MarkerSize', 8, 'DisplayName', 'Selected Peaks');
    end
    close all;

end