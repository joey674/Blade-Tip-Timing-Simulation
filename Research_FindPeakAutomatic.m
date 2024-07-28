function Research_FindPeakAutomatic(blade, EO)
    %{
        init params
    %}
    n_blades = length(blade);
                                                 
    %{
        find mode interval
        % 先计算所有叶片数据合并在一起;这里不归一化,虽然会有大数据罩住小数据的问题,但是如果归一化的话如果有一个纯噪音数据对整个模型的影响很大; 噪声会直接浮在整个归一化的图的上半部分
        % 再计算每一个freq点的最大值
        % 观测整体的模态 这个图会显示出所有可能的模态
        % 这里最后觉得先不高度平滑 不然有些模态会直接被抹除; 但是我还是想稍微平滑一点
    %}
    for i = 1:n_blades
        lens(i) = length([blade{i}.magn]);
    end
    magn_all = NaN(max(lens), n_blades);
    for i = 1:n_blades
        magn = smoothdata([blade{i}.magn], 'gaussian', 1000); 
        magn_all(1:length([blade{i}.magn]), i) = magn;
    end
    magn_max = max(magn_all, [], 2, "omitnan");
    % magn_max = smoothdata(magn_max, 'gaussian', 500); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{ 
        set parameter
        先用blade1的freq作为总共的freq 因为所有freq都一样
    %}
    freq = [blade{1}.freq]; 
    min_peak_prominence = 0.1 * (max(magn_max) - mean(magn_max)); 
    min_peak_height = 0.05 * (max(magn_max) - mean(magn_max)) + mean(magn_max); 
    min_peak_interval_width = 16;
    n_peaks = 10;   
    min_peak_prominence_tro = 0.01 * (max(magn_max) - mean(magn_max)); 

    %{
        find peaks 
        这里loc返回的就不是index了 是freq
    %}
    [pks, locs, widths, prominences] = findpeaks(magn_max, freq, ...
        'MinPeakProminence', min_peak_prominence, ...
        'MinPeakHeight', min_peak_height, ...
        'NPeaks', n_peaks...
    );

    %{
        find troughs
        这里loc_tro返回的就不是index了 是freq
    %}
    [tro_pks, tro_locs] = findpeaks(-magn_max, freq, ...
        'MinPeakProminence', min_peak_prominence_tro ...
    );

    %{
        calculate n_modes
        给后面使用的参数
    %}
    n_modes = length(pks);

    %{
        calculate peak interval
        每个区间的宽度将取 min_peak_interval_width 和波峰前后两个波谷距离之间的较小值
        首先就是模态之间不能有并集;
    %}
    peak_intervals = zeros(n_modes, 2);
    for i = 1:n_modes
        left_trough = max(tro_locs(tro_locs < locs(i)));
        right_trough = min(tro_locs(tro_locs > locs(i)));
        if isempty(left_trough) 
            left_trough = freq(1);
        end
        if isempty(right_trough)
            right_trough = freq(end);
        end

        if min_peak_interval_width/2 <= (right_trough - locs(i))
            peak_intervals(i, 2) = locs(i) + min_peak_interval_width/2;
        else
            peak_intervals(i, 2) = right_trough;
        end
        if min_peak_interval_width/2 <= (locs(i) - left_trough)
            peak_intervals(i, 1) = locs(i) - min_peak_interval_width/2;
        else
            peak_intervals(i, 1) = left_trough;
        end
    end

    %{
        plot
    %}
    fig = figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
    title(sprintf('EO%d', EO));
    xlim([13600 13900]);
    xlabel('Frequency ');
    ylabel('Magnitude ');
    set(gcf, 'WindowStyle');
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
             sprintf('F:%.2f\n H:%.2f\n P:%.2f\n W:%.2f',locs(i), pks(i), prominences(i), peak_intervals(i, 2) - peak_intervals(i, 1)), ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    legend('show');
    hold off;

    %% Save the figure
    saveas(fig, sprintf('mode_EO%d.png',EO));
    peak_freqs = locs;
    
    %{
        n_modes         所有叶片 总共有多少模态
        peak_intervals  所有模态的峰值处于的可能区间
    %}
    clearvars -except blade EO n_blades n_modes peak_intervals peak_freqs


    %% 对每个叶片执行寻找峰值策略
    % peaks_all(n_blades) = struct('peaks', struct('magn', {}, 'freq', {}, 'prominences', {}, 'widths', {}));
    % for blade_idx = 1:n_blades
    %     %{
    %         init and smooth
    %         平滑方式先用rlowess 之后修改成写好的自制方法
    %     %}
    %     blade_data = blade{blade_idx};
    %     freq = [blade_data.freq];
    %     magn = [blade_data.magn];
    %     phase = [blade_data.phase];
    %     err = [blade_data.err]; 
    %     magn = MDOF_ReduceNoise(magn, freq, err); 
    % 
    %     %% 根据所有叶片的相似模态形状寻找峰值 
    %     %{
    %         % set parameters
    %         这里由于没有过度平滑 所以这里的findpeak获得的参数prominence和width都没有意义 会被噪音干扰
    %     %}  
    %     peaks = struct('magn', {}, 'freq', {},'idx',{});
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %{
    %         % find peak
    %         % 对magn过度平滑, 寻找峰值
    %         % 判断峰值是否存在于模态区间内;如果存在则这个模态存在,记录这个模态区间
    %         % 对没有过度平滑的magn在此区间求最大值,记录freq
    %         % 更新 peaks 结构体数组，包含 magn, freq
    %     %}
    %     magn_oversmoothed = smoothdata(magn, 'gaussian', 1000); 
    %     min_peak_prominence = 0.03 * (max(magn_oversmoothed) - mean(magn_oversmoothed)); 
    %     min_peak_height = 0.06 * (max(magn_oversmoothed) - mean(magn_oversmoothed)) + mean(magn_oversmoothed); 
    %     [pks, locs, widths, prominences] = findpeaks(magn_oversmoothed, freq, ...
    %         'MinPeakProminence', min_peak_prominence, ...
    %         'MinPeakHeight', min_peak_height ...
    %     );
    % 
    % 
    %     %{
    %         % 现在开始遍历所有记录的模态区间,在未平滑的magn中找到该区间的最大值
    %         % 比较精确峰值与模态中的峰值,设置阈值判断是否接近如果找不到相近的模态
    %         峰值，则删除该峰值
    %     %}
    %     freq_diff_threshold = 2;
    %     for i = 1:n_modes
    %         interval = peak_intervals(i, :);
    %         disp(['Interval: [', num2str(interval(1)), ', ', num2str(interval(2)), ']']);
    % 
    %         in_interval = (locs >= interval(1)) & (locs <= interval(2));
    %         interval_locs = locs(in_interval);
    %         interval_pks = pks(in_interval);
    % 
    %         if isempty(interval_pks)
    %             disp('No peaks found in this interval.');
    %             continue;
    %         end
    % 
    %         in_interval_original = (freq >= interval(1)) & (freq <= interval(2));
    %         interval_freq = freq(in_interval_original);
    %         interval_magn = magn(in_interval_original);
    %         [max_magn, max_idx] = max(interval_magn);
    %         max_freq = interval_freq(max_idx);
    % 
    %         [min_diff, idx] = min(abs(peak_freqs(i) - max_freq));
    %         if min_diff < freq_diff_threshold 
    %             new_peak.magn = max_magn;
    %             new_peak.freq = max_freq;
    %             new_peak.idx = find(max_freq == freq);
    %             peaks = [peaks, new_peak];
    %         end 
    %     end
    % 
    %     %{
    %         plot
    %         % 对于每个blade绘制三个图:
    %         % subplot1是过度平滑的magn 并标注出峰值;
    %         % subplot2是原magn,并标注出记录的峰值
    %         % subplot3是最终确定的峰值
    %     %}
    %     figure('units', 'normalized', 'outerposition', [0 0 0.7 0.7]);
    %     set(gcf, 'WindowStyle', 'docked');
    % 
    % 
    %     subplot(2, 1, 1);
    %     plot(freq, magn_oversmoothed, 'Color', [0.7, 0.9, 1.0]);
    %     hold on;
    %     plot(locs, pks, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Detected Peak');
    %     for j = 1:length(locs)
    %         text(locs(j), pks(j), sprintf('P: %.2f\nW: %.2f', prominences(j), widths(j)), ...
    %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    %     end
    %     title(sprintf('EO%d blade%d - Oversmoothed with all found peaks', EO, blade_idx));
    %     xlabel('Frequency');
    %     ylabel('Magnitude');
    %     legend('show');
    %     hold off;
    % 
    %     subplot(2, 1, 2);
    %     plot(freq, magn, 'Color', [0.7, 0.9, 1.0], 'DisplayName', 'checked with mode plot');
    %     hold on;
    %     for j = 1:length(peaks)
    %         plot(peaks(j).freq, peaks(j).magn, 'ro', 'MarkerFaceColor', 'r');
    %         text(peaks(j).freq, peaks(j).magn, sprintf('F:%.2f\n M: %.2f', ...
    %             peaks(j).freq, peaks(j).magn), ...
    %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    %     end
    %     title(sprintf('EO%d blade%d - checked with mode plot', EO, blade_idx));
    %     xlabel('Frequency');
    %     ylabel('Magnitude');
    %     legend('show');
    %     hold off;
    % 
    %     %{ 
    %         获得结果数组
    %     %}
    %     peaks_all(blade_idx).peaks = peaks;
    % 
    %     %% 使用相位密度寻找峰值
    %     % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');
    %     % title(sprintf('EO%d, blade%d', EO, blade_idx));xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
    %     % hold on;
    %     % legend;  
    %     % 
    %     % subplot(4,1,1);
    %     % plot(freq, phase,'o','Color', [0.8, 0.9, 1.0]);  
    %     % 
    %     % % 定义频率范围和区间数
    %     % num_bins = 100;  %%%%%%%%%%%%%%%%%%%% 区间精度,分得越多得到的范围越精细
    %     % freq_edges = linspace(min(freq), max(freq), num_bins + 1);
    %     % phase_range = [-0.5, 0.5]; %%%%%%%%%%%%%%%%% 相位范围
    %     % 
    %     % % 初始化密度数组
    %     % phase_density = zeros(1, num_bins);
    %     % freq_centers = (freq_edges(1:end-1) + freq_edges(2:end)) / 2;
    %     % 
    %     % % 计算每个频率区间内相位在指定范围内的密度
    %     % for i = 1:num_bins
    %     %     % 当前频率区间
    %     %     idx = freq >= freq_edges(i) & freq < freq_edges(i+1);
    %     %     current_phase = phase(idx);
    %     % 
    %     %     % 计算相位在-0.5到0.5之间的点数比例作为密度
    %     %     num_points_in_range = sum(current_phase >= phase_range(1) & current_phase <= phase_range(2));
    %     %     total_points = length(current_phase);
    %     % 
    %     %     if total_points > 0
    %     %         phase_density(i) = num_points_in_range / total_points;
    %     %     else
    %     %         phase_density(i) = 0;
    %     %     end
    %     % end
    %     % 
    %     % % % 绘制相位密度图
    %     % subplot(4,1,2);
    %     % plot(freq_centers, phase_density, 'o', 'LineWidth', 2);
    %     % 
    %     % % 获取所有密度在threshold以下的频率区间的范围
    %     % low_density_threshold = 0.1;%%%%%%%%%%%%%%%%%%%%%%相位密度,中间空的地方就是有相变
    %     % low_density_ranges = [];
    %     % 
    %     % for i = 1:num_bins
    %     %     if phase_density(i) < low_density_threshold
    %     %         low_density_ranges = [low_density_ranges; freq_edges(i), freq_edges(i+1)];
    %     %     end
    %     % end
    %     % 
    %     % subplot(4,1,3);
    %     % 
    %     % % 降噪
    %     % magn = Damping_NoiseFilter(magn);
    %     % plot(freq, magn,'Color', [0.8, 0.9, 1.0]); hold on; 
    %     % 
    %     % % 
    %     % min_prominence = 1/4*mean(magn);%%%%%%%%%%%%%%%%%%%%%峰值的显著性 这个值可以小一点,多收入一些peak
    %     % min_peakdistance = 1/30*(max(freq)-min(freq));%%%%%%%%%%%%%%%%%%峰值之间的距离 这个值可以设置得大一些,避免peak挤在一起
    %     % [pks, locs] = findpeaks(magn, freq, 'MinPeakProminence', min_prominence,'MinPeakDistance',min_peakdistance);
    %     % % plot(locs, pks, 'o', 'DisplayName', 'Peaks through findpeak');
    %     % 
    %     % %
    %     % % 标注低密度频率区间(峰值所在区间)
    %     % for i = 1:size(low_density_ranges, 1)
    %     %     x_rect = [low_density_ranges(i, 1), low_density_ranges(i, 2), low_density_ranges(i, 2), low_density_ranges(i, 1)];
    %     %     y_rect = [min(magn), min(magn), max(magn), max(magn)];
    %     %     fill(x_rect, y_rect, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Low Density Region for peak');
    %     % end
    %     % 
    %     % % 筛选在低密度范围内的峰值
    %     % selected_pks = [];
    %     % selected_locs = [];
    %     % 
    %     % for i = 1:length(locs)
    %     %     for j = 1:size(low_density_ranges, 1)
    %     %         if locs(i) >= low_density_ranges(j, 1) && locs(i) < low_density_ranges(j, 2)
    %     %             selected_pks = [selected_pks, pks(i)];
    %     %             selected_locs = [selected_locs, locs(i)];
    %     %             break;
    %     %         end
    %     %     end
    %     % end
    %     % 
    %     % % 绘制筛选出的峰值
    %     % plot(selected_locs, selected_pks, 'bo', 'MarkerSize', 8, 'DisplayName', 'Selected Peaks');
    % end

end