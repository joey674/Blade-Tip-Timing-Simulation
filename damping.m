dataset = '../App4_SD24_37000_down_offen_220_run2_RPMRange-34200_34700_thres_S7_RotorMethod.mat';
n_blades = 12;
EO = 24;

%% choose the situation
load_data_form = "single_EO";
% load_data_form = "multi_EO";

% pre_process = "re_order_amplitude";
pre_process = "smooth_rotating_speed";

% method = "halfpower_bandwidth_method";
% method = "halfpower_bandwidth_method_third_correction";
method = "single_degree_of_freedom_approximation";
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[**********damping ratio calculation start.**********]\n');

%% load RPM and Magn from file
if (load_data_form == "single_EO")
    data = load(dataset, 'mean_RPM' , 'P_Magn');
    if isfield(data,'mean_RPM')
        excite_freq = data.mean_RPM * EO / 60;    
    else
        fprintf('Variable does not exist in the file\n');
    end
    if isfield(data, 'P_Magn')     
        Magn = cell(1, n_blades);     
        for i = 1:n_blades
            Magn{i} = data.P_Magn{i,24,:};  
        end    
    else
        fprintf('Variable does not exist in the file\n');
    end
end
if(load_data_form == "multi_EO")



end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pre process of data
if(pre_process == "re_order_amplitude")
    blade = cell(1, n_blades);   
    % find max length
    max_length = max([length(excite_freq), cellfun(@length, Magn)]);   
    for i = 1:n_blades
        % Assign Magn{i} and excite_freq to the structure
        blade{i} = repmat(struct('freq', NaN, 'magn', NaN), 1, max_length);       
        len_freq = length(excite_freq);
        len_magn = length(Magn{i});
        for j = 1:max_length
            if j <= len_freq
                blade{i}(j).freq = excite_freq(j);
            end
            if j <= len_magn
                blade{i}(j).magn = Magn{i}(j);
            end
        end    
        % Sort the structure array by freq
        [~, idx_sort] = sort([blade{i}.freq]);
        blade{i} = blade{i}(idx_sort);   
        % Remove elements contain nan
        idx_remove = isnan([blade{i}.freq]) | isnan([blade{i}.magn]);                   
        blade{i}(idx_remove) = [];   
    end
end
if(pre_process == "smooth_rotating_speed")
    blade = cell(1, n_blades);   
    for i = 1:n_blades
        % Remove NaNs from excite_freq
        valid_indices = ~isnan(excite_freq);
        magn_i = Magn{i}(valid_indices);
        freq_i = excite_freq(valid_indices);        
        % Find the shorter length
        len = min(length(magn_i), length(freq_i));        
        % Cut down magn_i and freq_i
        magn_i = magn_i(1:len);
        freq_i = freq_i(1:len);        
        % Sort freq_i
        [freq_i, idx_sort] = sort(freq_i);    
        % Assign magn_i and sorted freq_i to the structure
        blade{i} = struct('freq', num2cell(freq_i), 'magn', num2cell(magn_i));
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwidth_method)
if(method == "halfpower_bandwidth_method")    
    n_cols = 2;
    n_rows = ceil(n_blades / n_cols);
    num_blades = n_blades;
    colors = jet(num_blades);  
    for i = 1:num_blades
        figure('units','normalized','outerposition',[0 0 1 1]);
        % get magn and freq for each blade
        blade_data = blade{i};
        freq_i = [blade_data.freq];
        magn_i = [blade_data.magn]; 
        % plot origin diagram 
        yyaxis left;
        plot(freq_i, magn_i);
        ylabel('Magnitude');
        hold on;      
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess',1000);%%%%%%%%%%%%%%%%%     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10); %%%%%%%%%%%%%%%% 
        % plot smoothed magn
        plot(freq_i, Magn_each_blade, 'LineWidth', 2, 'Color', 'r');
        hold on;
        % deal each peak    
        damping_ratio = zeros(size(pks));
        for j = 1:length(pks)
            % find half bandwith point 
            f1 = freq_each_blade(find(Magn_each_blade(1:locs(j)) <= pks(j)/sqrt(2), 1, 'last'));
            f2 = freq_each_blade(find(Magn_each_blade(locs(j):end) <= pks(j)/sqrt(2), 1, 'first') + locs(j) - 1);
            %check if f1,f2 are found
            if isempty(f1)
                f1 = freq_each_blade(1);
            end
            if isempty(f2)
                f2 = freq_each_blade(end);
            end
            % calculate damping ratio
            damping_ratio(j) = (f2 - f1) / (2 * freq_each_blade(locs(j)));
            % plot peaks with different colors and mark peak index
            plot(freq_each_blade(locs(j)), pks(j), 'p', 'MarkerSize', 12, 'Color', colors(j, :));
            text(freq_each_blade(locs(j)), pks(j), ['Peak ', num2str(j)], 'Color', colors(j, :));
            % plot half power points
            plot(f1, pks(j)/sqrt(2), 'x', 'MarkerSize', 12, 'Color', 'g');
            plot(f2, pks(j)/sqrt(2), 'x', 'MarkerSize', 12, 'Color', 'g');          
        end
        % plot damping ratio
        yyaxis right;
        plot(freq_each_blade(locs), damping_ratio, 'o', 'MarkerSize', 12);
        % Add damping ratio and corresponding freq value with loop
        for j = 1:length(damping_ratio)
            text(freq_each_blade(locs(j)), damping_ratio(j), ['DR: ', num2str(damping_ratio(j)), newline, 'Freq: ', num2str(freq_each_blade(locs(j)))], 'Color', colors(j, :));
        end
        ylabel('Damping Ratio');
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
        title(['Blade ', num2str(i)]);
        hold off;
        % save file
        saveas(gcf, fullfile('..', 'DampingRatio', ['Blade_', num2str(i), '_HPBWM','.png']));
        close(gcf); 
    end        
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method_third_correction)
if(method == "halfpower_bandwidth_method_third_correction")
    n_cols = 2;
    n_rows = ceil(n_blades / n_cols);
    num_blades = n_blades; 
    colors = jet(num_blades);  
    for i = 1:num_blades  
        figure('units','normalized','outerposition',[0 0 1 1]);      
        % get magn and freq for each blade
        blade_data = blade{i};
        freq_i = [blade_data.freq];
        magn_i = [blade_data.magn]; 
        % plot origin diagram 
        yyaxis left;
        plot(freq_i, magn_i);
        ylabel('Accelaration');
        hold on;        
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        % Convert to acceleration response function
        Magn_each_blade = (freq_each_blade).^2 .* Magn_each_blade;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess',1000); %%%%%%%%%%%%%%    
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',1500000000);%%%%%%%%%%%%%   
        % plot smoothed magn
        plot(freq_i, Magn_each_blade, 'LineWidth', 2, 'Color', 'r');
        hold on;
        % deal each peak    
        damping_ratio = zeros(size(pks));
        for j = 1:length(pks)
            % find half bandwidth point 
            f1 = freq_each_blade(find(Magn_each_blade(1:locs(j)) <= pks(j)/sqrt(2), 1, 'last'));
            f2 = freq_each_blade(find(Magn_each_blade(locs(j):end) <= pks(j)/sqrt(2), 1, 'first') + locs(j) - 1);
            %check if f1,f2 are found
            if isempty(f1)
                f1 = freq_each_blade(1);
            end
            if isempty(f2)
                f2 = freq_each_blade(end);
            end
            % calculate damping ratio
            sol = roots([8,0,2,(f1 - f2) / freq_each_blade(locs(j))]); 
            damping_ratio(j) = sol(imag(sol) == 0 & sol >= 0 & sol <= 1);
            % plot peaks with different colors and mark peak index
            plot(freq_each_blade(locs(j)), pks(j), 'p', 'MarkerSize', 12, 'Color', colors(j, :));
            text(freq_each_blade(locs(j)), pks(j), ['Peak ', num2str(j)], 'Color', colors(j, :));
            % plot half power points
            plot(f1, pks(j)/sqrt(2), 'x', 'MarkerSize', 12, 'Color', 'g');
            plot(f2, pks(j)/sqrt(2), 'x', 'MarkerSize', 12, 'Color', 'g');
        end
        % plot damping ratio
        yyaxis right;
        plot(freq_each_blade(locs), damping_ratio, 'o', 'MarkerSize', 12);
        % Add damping ratio and corresponding freq value with loop
        for j = 1:length(damping_ratio)
            text(freq_each_blade(locs(j)), damping_ratio(j), ['DR: ', num2str(damping_ratio(j)), newline, 'Freq: ', num2str(freq_each_blade(locs(j)))], 'Color', colors(j, :));
        end
        ylabel('Damping Ratio');
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
        title(['Blade ', num2str(i)]);
         hold off;
        % save file
        saveas(gcf, fullfile('..', 'DampingRatio', ['Blade_', num2str(i),'_HPBWM3C', '.png']));
        close(gcf); 
    end        
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(single_degree_of_freedom_approximation)
if(method == "single_degree_of_freedom_approximation")    
    n_cols = 2;
    n_rows = ceil(n_blades / n_cols);
    num_blades = n_blades;
    colors = jet(num_blades);  
    for i = 1:num_blades
        figure('units','normalized','outerposition',[0 0 1 1]);
        % get magn and freq for each blade/plot origin 
        blade_data = blade{i};
        freq_i = [blade_data.freq];
        magn_i = [blade_data.magn]; 
        yyaxis left;
        plot(freq_i, magn_i);
        ylabel('Magnitude');
        hold on;   
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;

        %smooth data to find peak and trough 
        Magn_smoothed = smoothdata(Magn_each_blade,'loess',1000);
        [peaks,peaks_locs] = findpeaks(Magn_smoothed,'MinPeakProminence',10); 
        [troughs, trough_locs] = findpeaks(-Magn_smoothed,'MinPeakProminence',10);
        troughs = -troughs;  

        % 设定满足条件的拟合差值 
        satisfied_error = 1;
        
        % 初始化空数组来存储每个波峰的最佳数据范围和误差
        best_ranges = zeros(length(peaks_locs), 2);
        best_errors = inf(length(peaks_locs), 1);
        
        % 处理每个波峰 会把得到的最佳数据集和拟合最好的误差返回出来.
        % 但我这里其实还要返回阻尼比
        for i = 1:length(peaks_locs)
            [best_ranges(i, :), best_errors(i)] = fitSDOF(peaks_locs, trough_locs, i, Magn_each_blade, satisfied_error);
        end






    end        
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[**********damping ratio calculation finish.*********]\n');

%% used funcitons
function H = sdof_response(w, w_n, xi)
    H = 1 ./ sqrt((1 - (w / w_n).^2).^2 + (2 * xi * (w / w_n)).^2);
end
function L = loss_function(p, w, Magn)
    % p(1) is natural frequency, p(2) is damping ratio
    H_model = sdof_response(w, p(1), p(2));
    L = sum((Magn - H_model).^2);
end
function [best_range, best_error] = fitSDOF(peaks_locs, trough_locs, i, Magn_each_blade, satisfied_error)
    % 初始化最佳范围和误差
    best_range = [];
    best_error = inf;    
    % 初始化数据集范围 先从最近的左右波谷开始
    left_trough = max(trough_locs(trough_locs < peaks_locs(i)));
    right_trough = min(trough_locs(trough_locs > peaks_locs(i)));    
    while true
        % 使用当前的数据范围进行拟合
        current_range = [left_trough, right_trough];   

        %todo
        current_error = 1; 
        

        % 如果当前误差小于best_error，更新best_error和best_range
        if current_error < best_error
            best_error = current_error;
            best_range = current_range;
        end
         % 如果达到满足误差或遍历完所有波谷，跳出循环
        if current_error <= satisfied_error || isempty(left_trough) && isempty(right_trough)
            break;
        end
        % 扩大数据范围:选择当前谷的下一个谷,但是不会越过最近的波峰
        left_trough = max(trough_locs(trough_locs < left_trough & (i == 1 || trough_locs > peaks_locs(i-1))));
        right_trough = min(trough_locs(trough_locs > right_trough & (i == length(peaks_locs) || trough_locs < peaks_locs(i+1))));
    end
end    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




