% include using function damping_re_smooth(Magn_each_blade,locs)

dataset = '../App4_SD24_37000_down_offen_220_run2_RPMRange-34200_34700_thres_S7_RotorMethod.mat';
n_blades = 12;
EO = 24;

%% choose a method or plot struct
method = "default_as_non";
% method = "halfpower_bandwith_method";
% method = "halfpower_bandwith_method_v1";
method = "halfpower_bandwith_method_v2";
%method = "halfpower_bandwith_method_third_correction";

plot_struct = "default_as_all";
% plot_struct = "debug_mode_show_first";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[**********damping ratio calculation start.**********]\n');

%% load RPM and Magn from file
data = load(dataset, 'mean_RPM' , 'P_Magn');
if isfield(data,'mean_RPM')
    excite_freq = data.mean_RPM * EO;    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method)
if(method == "halfpower_bandwith_method")
    idx_remove = isnan(excite_freq);
    excite_freq(idx_remove) = [];
    for i = 1:n_blades
        Magn{i}(idx_remove) = [];
    end
    len = length(excite_freq);
    for i = 1:n_blades
        Magn{i} = Magn{i}(1:len);
    end
    % sort excite_freq 
    [excite_freq, idx_sort] = sort(excite_freq);
    % reorder Magn 
    for i = 1:n_blades
        Magn{i} = Magn{i}(idx_sort);
    end

    damping_ratios = cell(1, n_blades);
    for i = 1:n_blades     
        Magn_each_blade = Magn{i};
        freq_each_blade = excite_freq;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess');     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10);  
        %re-smooth data
        %Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
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
        end
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method_v1)
if(method == "halfpower_bandwith_method_v1")
    % pre process of data
    idx_remove = isnan(excite_freq);
    excite_freq(idx_remove) = [];
    for i = 1:n_blades
        Magn{i}(idx_remove) = [];
    end
    len = length(excite_freq);
    for i = 1:n_blades
        Magn{i} = Magn{i}(1:len);
    end
    % sort excite_freq 
    [excite_freq, idx_sort] = sort(excite_freq);
    % reorder Magn 
    for i = 1:n_blades
        Magn{i} = Magn{i}(idx_sort);
    end
    
    %calculate damping ratios
    figure('units','normalized','outerposition',[0 0 1 1]);
    % n_cols = 2;
    % n_rows = ceil(n_blades / n_cols);
    % for i = 1:n_blades
        n_cols = 1;n_rows = 1;for i = 1:1%(one blade plot uncomment this)    
        subplot(n_rows, n_cols, i);
        % get magn and freq for each blade
        magn_i = Magn{i};
        freq_i = excite_freq;    
        % plot origin diagram 
        yyaxis left;
        plot(freq_i, magn_i);
        ylabel('Magnitude');
        hold on;        
        %%%%%%%%%%%%%%%%%%%%%
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess');     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10);  
        %re-smooth data
        % Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
        % plot smoothed magn
        plot(freq_i, Magn_each_blade, 'LineWidth', 2, 'Color', 'r');
        hold on;
        %%%%%%%%%%%%%%%%%%%%
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
            % plot peaks
            plot([freq_each_blade(locs(j)), freq_each_blade(locs(j))], [min(magn_i), max(magn_i)], 'r-');
        end
        % plot damping ratio
        yyaxis right;
        plot(freq_each_blade(locs), damping_ratio, 'o');
        ylabel('Damping Ratio');
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
    end        
    title(['Blade ', num2str(i)]);
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method_v2)
if(method == "halfpower_bandwith_method_v2")    
    % pre process of data
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
    %%%%%%%%%%%%%%%%%%%%%

    %calculate damping ratios
    figure('units','normalized','outerposition',[0 0 1 1]);
    if (plot_struct == "default_as_all")    
        n_cols = 2;
        n_rows = ceil(n_blades / n_cols);
        num_blades = n_blades;
    end
    if (plot_struct == "debug_mode_show_first") 
         n_cols = 1;
         n_rows = 1;
         num_blades = 1;
    end    
    for i = 1:num_blades  
        subplot(n_rows, n_cols, i);
        % get magn and freq for each blade
        blade_data = blade{i};
        freq_i = [blade_data.freq];
        magn_i = [blade_data.magn]; 
        % plot origin diagram 
        yyaxis left;
        plot(freq_i, magn_i);
        ylabel('Magnitude');
        hold on;        
        %%%%%%%%%%%%%%%%%%%%%
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess');     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10);  
        %re-smooth data
        % Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
        % plot smoothed magn
        plot(freq_i, Magn_each_blade, 'LineWidth', 2, 'Color', 'r');
        hold on;
        %%%%%%%%%%%%%%%%%%%%
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
            % plot peaks
            plot([freq_each_blade(locs(j)), freq_each_blade(locs(j))], [min(magn_i), max(magn_i)], 'r-');
        end
        % plot damping ratio
        yyaxis right;
        plot(freq_each_blade(locs), damping_ratio, 'o');
        ylabel('Damping Ratio');
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
        title(['Blade ', num2str(i)]);
    end        
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method_third_correction)
if(method == "halfpower_bandwith_method_third_correction")
    damping_ratios = cell(1, n_blades);
    locs_ratios = cell(1, n_blades);  % Store locations for each blade's damping ratios
    colors = lines(n_blades); % generate a matrix of different colors
    figure;
    subplot(2,1,1); % first subplot for acceleration response
    hold on;
    for i = 1:n_blades
        Magn_each_blade = Magn{i};
        freq_each_blade = excite_freq;
        % Convert to acceleration response
        Magn_each_blade = (2 * pi * freq_each_blade).^2 .* Magn_each_blade;
        % Smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess');     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10);  
        %re-smooth data
        Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
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
            damping_ratio(j) = ((f2 - f1) / freq_each_blade(locs(j))) / (2 + 8 * (f2 - f1 / freq_each_blade(locs(j)))^2);
        end
        % store calculated damping ratio of all peaks for one blade 
        damping_ratios{i} = damping_ratio;
        locs_ratios{i} = locs; % Store corresponding locations

        % plot acceleration response for each blade
        x = 1:length(Magn_each_blade);
        plot(x*EO, Magn_each_blade, 'LineWidth', 2, 'Color', colors(i,:));
    end
    xlabel('Frequency');
    ylabel('Acceleration Response');
    title('Acceleration Response for All Blades');
    hold off;

    % plot damping ratios for all blades
    subplot(2,1,2); % second subplot for damping ratios
    hold on;
    for i = 1:n_blades
        scatter(locs_ratios{i}*EO, damping_ratios{i}, 20, colors(i,:), 'filled');  % plot with scatter
    end
    xlabel('Frequency');
    ylabel('Damping Ratio');
    title('Damping Ratios for All Blades');
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf('[**********damping ratio calculation finish.*********]\n');




