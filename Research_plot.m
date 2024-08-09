%% VA对比
function Research_plot(blade, EO)
    % 加载fig文件
    fig1 = openfig('Graph/EO24_MDOF_blade2_Set1.fig', 'invisible');
    fig2 = openfig('Graph/EO24_MDOF_blade2_Set2.fig', 'invisible');
    fig3 = openfig('Graph/EO24_MDOF_blade2_Set3.fig', 'invisible');
    fig4 = openfig('Graph/EO24_MDOF_blade2_Set4.fig', 'invisible');
    fig5 = openfig('Graph/EO24_MDOF_blade2_Set5.fig', 'invisible');

    % 创建一个新的图形窗口
    set(0,'DefaultFigureWindowStyle','normal');
    fig = figure('NumberTitle', 'off', 'Position', [50, 50, 600, 1200]);
    
    % 使用tiledlayout来创建一个5行1列的布局
    t = tiledlayout(5, 1, 'TileSpacing', 'none', 'Padding', 'none');

    % 对每个子图执行copyobj
    figs = {fig1, fig2, fig3, fig4, fig5}; % 将fig对象放入一个cell array中
    for blade_idx = 1:5
        ax = nexttile; % 创建下一个子图
        % 获取当前fig的轴对象
        current_ax = gca(figs{blade_idx});
        % 复制当前轴对象的内容到新的子图
        copyobj(allchild(current_ax), ax);

        % 格式化每个子图
        grid on;
        hold on;
        xlim([13680,13880]);
        ylim([0, 1]);

        % 设置纵坐标标题但隐藏刻度值
        ylabel(sprintf('Set %d ', blade_idx), "Rotation", 0);
        yticks([]);
        
        % 隐藏除最后一个子图以外的所有子图的x轴
        if blade_idx < 5
            ax.XAxis.Visible = 'off'; % 隐藏x轴刻度
        else
            xlabel("Frequency (Hz)"); % 设置最后一个子图的x轴标签
        end
    end

    % 添加图例
    % lgd = legend(ax, 'Amplitude', 'Error', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
    lgd.Layout.Tile = 'North'; % 将图例放在图形顶部

    % 输出EO信息
    disp(['EO', num2str(EO)]);
end



%% FindPeak StackPlot
% function Research_plot(blade, EO)
%     n_blades = length(blade);
% 
%     set(0,'DefaultFigureWindowStyle','normal');
%     fig = figure( 'NumberTitle', 'off', 'Position', [50, 50, 600, 1200]);
%     t = tiledlayout(fig, n_blades, 1);
% 
%     if EO == 24
%         x_left = 13680;
%         x_right = 13880;
%         normalized_factor = 120;
%     elseif EO == 20
%          x_left = 9420;
%         x_right = 9730;
%         normalized_factor = 55;
%     elseif EO == 8
%          x_left = 3770;
%         x_right = 3890;
%         normalized_factor = 65;
%     end
% 
%     peaks_idx_all = FindPeakAutomatic(blade);
% 
%     for blade_idx = n_blades:-1:1
%         % 初始化
%         fprintf('blade:%d\n', blade_idx); 
%         blade_data = blade{blade_idx};
%         freq = [blade_data.freq];
%         magn = [blade_data.magn];
%         phase = [blade_data.phase];
%         err = [blade_data.err]; 
% 
%         % 归一化
%         err = err / normalized_factor;
%         magn = magn / normalized_factor;
% 
%         % 降噪和平滑处理
%         magn = ReduceNoise(magn, freq, err);
%         err = smoothdata(err, "movmean", 200);
% 
%         % 在网格中绘制数据
%         ax(blade_idx) = nexttile;
%         grid on;
%         hold on;
% xlim([x_left,x_right]);
%         ylim([0, 1]);
% 
%         % 设置纵坐标标题但隐藏刻度值
%         ylabel(sprintf('Blade %d ', blade_idx), "Rotation", 0);
%         yticks([]);
% 
%         % 绘制幅值和误差
%         plot(freq, magn, 'LineWidth', 1.5, 'DisplayName', sprintf('Blade %d Amplitude', blade_idx));
%         plot(freq, err, 'LineWidth', 1.5, 'DisplayName', sprintf('Blade %d Error', blade_idx));
% 
%         %绘制峰值       
%         peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
%         % [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx);
%         peaks_idx = sort(peaks_idx);
%         if isempty(peaks_idx) 
%             disp('ERROR: can not find peaks for this blade, will skip');
%             continue;
%         end
%         plot(freq(peaks_idx), magn(peaks_idx), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%         hold off;
%     end
% 
%     % 添加图例和标签
%     lgd = legend(ax(n_blades), 'Amplitude', 'Error', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
%     lgd.Layout.Tile = 'North';
%     xlabel(ax(1), "Frequency (Hz)");
% 
%     % 隐藏除第一个以外的所有子图的 x 轴
%     for i_axis = 2:n_blades
%         ax(i_axis).XAxis.Visible = 'off';
%     end
% 
%     % 图形格式
%     t.TileSpacing = 'none';
% 
%     disp(['EO',num2str(EO)]);
% end

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


