function Research_Mode(blade,EO)
    %% init params
    n_blades = length(blade);

        figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle');
        hold on;
        legend;  

    %% deal every blades
    for blade_idx = 1:n_blades
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 

        freq = freq(:);

        % 找到频率数据的范围
        min_freq = floor(min(freq));
        max_freq = ceil(max(freq));
    
        % 初始化每1Hz区间的数据点计数
        bin_counts = zeros(max_freq - min_freq, 1);
    
        % 统计每1Hz区间的数据点数
        for i = min_freq:max_freq-1
            bin_counts(i - min_freq + 1) = sum(freq >= i & freq < i + 1);
        end
    
        % 创建1Hz区间标签
        bins = min_freq:max_freq-1;
    
        % 可视化每1Hz区间的数据点数
        figure;
        bar(bins, bin_counts, 'FaceColor', [0.2 0.2 0.5], 'EdgeColor', [0.2 0.2 0.5]);
        xlabel('Frequency (Hz)');
        ylabel('Number of Data Points');
        title('Number of Data Points per 1Hz Frequency Interval');
        grid on;

       
    end
end