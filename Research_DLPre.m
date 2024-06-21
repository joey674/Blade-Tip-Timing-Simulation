clear;
datasets = {
    '../Dataset/EO24/App4_SD24_37000_down_offen_220_run3_RPMRange-34100_34700_thres_S7_RotorMethod_movmean_10.mat', 12, [24];
    % '../Dataset/EO08_20/App4_DRH30_29300_down_offen_180_run1_RPMRange-28200_29200_thres_S8_RotorMethod_movmean_10.mat', 12, [20, 8];
    % '../Dataset/EO16/App4_SD24_37000_down_offen_220_run3_RPMRange-35600_36200_thres_S7_RotorMethod_movmean_10.mat', 12, [16];
    % '../Dataset/EO17/App4_SD24_37000_down_offen_220_run3_RPMRange-33600_34200_thres_S7_RotorMethod_movmean_10.mat', 12, [17];
    % '../Dataset/EO30/App4_DRH30_29300_down_offen_180_run1_RPMRange-27200_27700_thres_S8_RotorMethod_movmean_10.mat', 12, [30];
};

P = [];
T = [];

for k = 1:size(datasets, 1)
    dataset_path = datasets{k, 1};
    n_blades = datasets{k, 2};
    EO = datasets{k, 3};
    data = load(dataset_path, 'mean_RPM', 'P_Magn','P_Phase', 'Fit_Error');
    blade_sets = Damping_PreProcess(EO, data, n_blades);
    for i = 1:length(blade_sets)
        n_blades = length(blade_sets{i}); 
        target_length = 1000;
    
        for blade_idx = 1:n_blades
            disp(['blade:', num2str(blade_idx)]); % 载入数据
            blade_data = blade_sets{i}{blade_idx};
            freq = [blade_data.freq];
            magn = [blade_data.magn];
            phase = [blade_data.phase];
            err = [blade_data.err];
            magn = MDOF_ReduceNoise(magn, freq, err); % 降噪

            % 重采样 将数据降采样到目标长度
            original_length = length(magn);
            new_indices = linspace(1, original_length, target_length);
            magn = interp1(1:original_length, magn, new_indices, 'linear');
            freq = interp1(1:original_length, freq, new_indices, 'linear');
            phase = interp1(1:original_length, phase, new_indices, 'linear');
            err = interp1(1:original_length, err, new_indices, 'linear');
            P = [P; magn];

            
            peak_vector = zeros(1, length(magn));% 初始化T向量
            fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);% 图形设置
            set(gcf, 'WindowStyle', 'docked');
            title(sprintf('EO%d, blade%d', EO(1), blade_idx));
            xlabel('Index');
            ylabel('Magnitude');
            hold on;
            legend;
            plot(magn, 'Color', [0.8, 0.9, 1.0]);
            while true
                [x, ~] = ginput(1);
                if isempty(x)
                    break;
                end
                x = round(x);
                plot(x, magn(x), 'ro'); % 标注红点
                [x2, ~] = ginput(1);
                if isempty(x2)
                    break;
                end
                x2 = round(x2);
                plot(x2, magn(x2), 'bo'); % 标注蓝点
                % 确保 x < x2
                if x > x2
                    [x, x2] = deal(x2, x);
                end
                peak_vector(x:x2) = 1; % 标记区间内的T为1
            end
            T = [T; peak_vector];
            close(fig); % 关闭当前图形窗口
        end
    end
end

save('DL/peak.mat', 'P', 'T');




