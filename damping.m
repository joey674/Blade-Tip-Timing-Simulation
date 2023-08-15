dataset = '../App4_SD24_37000_down_offen_220_run2_RPMRange-34200_34700_thres_S7_RotorMethod.mat';
n_blades = 12;
EO = 24;

%% choose the situation
load_data_form = "single_EO";
% load_data_form = "multi_EO";

pre_process = "re_order_amplitude";
% pre_process = "smooth_rotating_speed";

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
    num_blades = n_blades;

    % deal with each blade
    for blade_idx = 1:num_blades
        figure('units','normalized','outerposition',[0 0 0.7 0.7]);   
        hold on;

        % get blade data
        blade_data = blade{blade_idx};
        freq_each_blade = [blade_data.freq];
        Magn_each_blade = [blade_data.magn]; 

        % plot origin graph
        plot(freq_each_blade, Magn_each_blade, 'k', 'DisplayName', 'Original Data');

        % find peaks and troughs
        Magn_smoothed = smoothdata(Magn_each_blade,'loess',100);
        [peaks,peaks_locs] = findpeaks(Magn_smoothed,'MinPeakProminence',10); 
        [troughs, trough_locs] = findpeaks(-Magn_smoothed,'MinPeakProminence',1);
        troughs = -troughs;  
        % plot peaks and troughs
        plot(freq_each_blade(peaks_locs), peaks, 'ro', 'DisplayName', 'Peaks');
        plot(freq_each_blade(trough_locs), troughs, 'bo', 'DisplayName', 'Troughs');
        legend;
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');

        % set satisfied approximate degree(let it be zero so never satisfied)
        satisfied_error = 0;     

        % init vector to store value
        best_ranges = zeros(length(peaks_locs), 2);
        best_errors = inf(length(peaks_locs), 1);    

        % deal every peaks
        for i = 1:length(peaks_locs)
            % init dataset to store best result
            best_range = [];
            best_error = inf;    

            % choose first pair of troughs
            left_trough = max(trough_locs(trough_locs < peaks_locs(i)));
            right_trough = min(trough_locs(trough_locs > peaks_locs(i)));  
            iteration = 0;
            
            % deal every pair of troughs
            while true                                             
                % init
                current_range = [left_trough, right_trough];    
                current_magn = Magn_smoothed(current_range(1):current_range(2));
                % current_magn = Magn_each_blade(current_range(1):current_range(2));
                current_freq = freq_each_blade(current_range(1):current_range(2));
                w_n = freq_each_blade(peaks_locs(i));
                w = current_freq;

                % set value for fmincon
                down_boundary = [0,0]; 
                upper_boundary = [inf,0.1];
                startvalue = [1/peaks(i),0.001];% 一个是激励因子的初始值,一个是阻尼比的初始值

                % 开始拟合 得到这个范围下最佳的阻尼比 以及这个阻尼比对应的error(error
                % 衡量这个最佳的阻尼比是否满足我们的需求)
                % start to aproximate
                options = optimoptions('fmincon', 'Display', 'iter');
                optimal_params = fmincon(@(params) computeError(params, w, w_n, current_magn), ...
                startvalue, [], [], [], [], down_boundary, upper_boundary, [], options);
                A = optimal_params(1);
                xi = optimal_params(2);
                

                % use DR that we get to calculate error              
                X = A ./ sqrt((1 - (w / w_n) .^ 2) .^ 2 + (2 * xi * (w / w_n)) .^ 2);
                current_error = computeError(optimal_params, w, w_n, current_magn);

                % check if satisfied then break
                if current_error < best_error
                    best_error = current_error;
                    best_range = current_range;   
                    % plot graph
                    fprintf("================\n");
                    color_depth = max(0, 1 - iteration*0.2);
                    plot(w, X, 'Color', [0, color_depth, 0], 'LineWidth',...
                    1+iteration*0.5);
                    fprintf("++++++++++++++++\n");
                    if best_error <= satisfied_error 
                        break;
                    end
                end  
                iteration = iteration + 1;

                % iterate troughs
                if i == 1 
                    valid_left_troughs = trough_locs(trough_locs < left_trough);
                else 
                    valid_left_troughs = trough_locs(trough_locs < left_trough & trough_locs > peaks_locs(i-1));
                end
                left_trough = max(valid_left_troughs);               
                if i == length(peaks_locs) 
                    valid_right_troughs = trough_locs(trough_locs > right_trough);
                else 
                    valid_right_troughs = trough_locs(trough_locs > right_trough & trough_locs < peaks_locs(i+1));
                end
                right_trough = min(valid_right_troughs);
                % its the last trough,break
                if isempty(left_trough) || isempty(right_trough)
                    break;
                end

            end
            % store result 
            best_ranges(i, :) = best_range;
            best_errors(i) = best_error;

        end
        
        % add plot info
        legend;
        title(['Blade ', num2str(blade_idx)]);
        xlabel('Frequency');
        ylabel('Magnitude');
        hold off;

    end        
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[**********damping ratio calculation finish.*********]\n');

%% funcitons for SDOF
function error = computeError(params, w, w_n, Magn_each_blade)
    A = params(1);
    xi = params(2);
    X = A ./ sqrt((1 - (w / w_n) .^ 2) .^ 2 + (2 * xi * (w / w_n)) .^ 2);
    len =length(w);
    % error = sum((Magn_each_blade - X) .^ 2)/len;
    error = norm(Magn_each_blade - X)/sqrt(len);
end 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




