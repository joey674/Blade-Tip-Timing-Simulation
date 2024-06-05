function Research_Mode(blade,EO)
    %% init params
    n_blades = length(blade);

    %% deal every blades
    for blade_idx = 1:n_blades
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 
      
        %%
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');
        title(sprintf('EO%d, blade%d', EO, blade_idx));
        hold on;
        legend;  

        % subplot(4,1,1);
        % plot(freq, magn,'o','Color', [0.8, 0.9, 1.0]);  hold on;
        % plot(freq, smoothdata(magn,'rlowess',200),'o');  
        % 
        % subplot(4,1,2);
        % plot(freq, phase,'o','Color', [0.8, 0.9, 1.0]);  
        
        phase_smoothed = smoothdata(phase,'movmean',200);
        subplot(2,1,1);
        plot(freq, phase_smoothed,'o','Color', [0.8, 0.9, 1.0]);  

        dphase_dfreq = diff(phase_smoothed);
        subplot(2,1,2);
        plot(dphase_dfreq, 'r-', 'DisplayName', 'dPhase/dFreq');


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
end