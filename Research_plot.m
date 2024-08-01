%% stackplot
function Research_plot(blades, EO)
    n_blades = length(blades);

    set(0,'DefaultFigureWindowStyle','normal');
    titel = sprintf('EO_%d_Stacked_Plot', EO);
    fig = figure('Name', titel, 'NumberTitle', 'off', 'Position', [50, 50, 1600, 1200]);
    t = tiledlayout(fig, n_blades, 1);

    for blade_idx = n_blades:-1:1
        % 初始化
        fprintf('blade:%d\n', blade_idx); 
        blade_data = blades{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 

        % 归一化
        % magn = magn / max(magn);
        % err = err / max(magn);

        % 降噪和平滑处理
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, "movmean", 200);

        % 在网格中绘制数据
        ax(blade_idx) = nexttile;
        grid on;
        hold on;
        xlim([13680, 13880]);
        % ylim([0, 1]);
        
        % 设置纵坐标标题但隐藏刻度值
        ylabel(sprintf('Blade %d ', blade_idx), "Rotation", 0);
        yticks([]);

        % 绘制幅值和误差
        plot(freq, magn, 'LineWidth', 1.5, 'DisplayName', sprintf('Blade %d Amplitude', blade_idx));
        plot(freq, err, 'LineWidth', 1.5, 'DisplayName', sprintf('Blade %d Error', blade_idx));
        hold off;
    end

    % 添加图例和标签
    lgd = legend(ax(n_blades), 'Amplitude', 'Error', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
    lgd.Layout.Tile = 'North';
    xlabel(ax(1), "Frequency (Hz)");

    % 隐藏除第一个以外的所有子图的 x 轴
    for i_axis = 2:n_blades
        ax(i_axis).XAxis.Visible = 'off';
    end

    % 图形格式
    t.TileSpacing = 'none';
end



%% 降噪
% function Research_plot(blades, EO)
%     %% init params
%     n_blades = length(blades);
% 
%     %% get peaks_idx    
%     peaks_idx_all = FindPeakAutomatic(blades);
% 
%     %% deal with every blade
%     for blade_idx = 1:n_blades
%         %% init
%         fprintf('blade:%d\n', blade_idx); 
%         blade_data = blades{blade_idx};
%         freq = [blade_data.freq];
%         magn = [blade_data.magn];
%         phase = [blade_data.phase];
%         err = [blade_data.err]; 
%         % normalized
%         magn = magn/max(magn);
%         err = err/max(magn);
% 
%         %% plot the original data
%         figure;
%         set(gcf, 'WindowStyle', 'docked');
%         hold on;
% 
%         %% apply noise reduction methods
%         windowSize = 400;
% 
%         % movmean
%         magn_movmean = movmean(magn, windowSize);
%         plot(freq, magn_movmean, 'DisplayName', 'Movmean', 'LineWidth', 1.5);
% 
%         % % movmedian
%         % magn_rloess = smoothdata(magn, 'movmedian', windowSize);
%         % plot(freq, magn_rloess, 'DisplayName', 'movmedian', 'LineWidth', 1.5);
% 
%          % rloess
%         magn_rloess = smoothdata(magn, 'rloess', windowSize);
%         plot(freq, magn_rloess, 'DisplayName', 'rloess', 'LineWidth', 1.5);
% 
%         % sgolay
%         magn_sgolay = smoothdata(magn, 'sgolay', windowSize);
%         plot(freq, magn_sgolay, 'DisplayName', 'Sgolay', 'LineWidth', 1.5);
% 
%         %% configure the plot
%         xlim([13680, 13880]);
%         ylim([0, 1]);
%         xlabel('Frequency (Hz)');
%         ylabel('Normalized Amplitude');
%         legend('show');
%         hold off;
%     end
% end
% 


