dataset = '../App4_SD24_37000_down_offen_220_run2_RPMRange-34200_34700_thres_S7_RotorMethod.mat';
n_blades = 12;
EO = 24;

%% choose the situation
load_data_form = "single_EO";
% load_data_form = "multi_EO";

pre_process = "re_order_amplitude";
% pre_process = "smooth_rotating_speed";

method = "non";
method = "halfpower_bandwidth_method";
% method = "halfpower_bandwidth_method_third_correction";
% method = "single_degree_of_freedom_approximation";

plot_struct = "one";
plot_struct = "four";
% plot_struct = "all";

data_output = "together"
% data_output = "seperate"
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


end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwidth_method)
if(method == "halfpower_bandwidth_method")    
    figure('units','normalized','outerposition',[0 0 1 1]);
    if (plot_struct == "all")    
        n_cols = 2;
        n_rows = ceil(n_blades / n_cols);
        num_blades = n_blades;
    end
    if (plot_struct == "one") 
         n_cols = 1;
         n_rows = 1;
         num_blades = 1;
    end    
    if (plot_struct == "four") 
         n_cols = 2;
         n_rows = ceil(4 / n_cols);
         num_blades = 4;
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
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess',1000);     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',10);  
        %re-smooth data
        % Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(halfpower_bandwith_method_third_correction)
if(method == "halfpower_bandwidth_method_third_correction")
    figure('units','normalized','outerposition',[0 0 1 1]);
    if (plot_struct == "all")    
        n_cols = 2;
        n_rows = ceil(n_blades / n_cols);
        num_blades = n_blades;
    end
    if (plot_struct == "one") 
         n_cols = 1;
         n_rows = 1;
         num_blades = 1;
    end    
    if (plot_struct == "four") 
         n_cols = 2;
         n_rows = ceil(4 / n_cols);
         num_blades = 4;
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
        ylabel('Accelaration');
        hold on;        
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        % Convert to acceleration response function%%%%%%%%%%%%%%%%%%%%%%
        Magn_each_blade = (freq_each_blade).^2 .* Magn_each_blade;
        %smooth data
        Magn_each_blade = smoothdata(Magn_each_blade,'loess',100);     
        % find peaks
        [pks,locs] = findpeaks(Magn_each_blade,'MinPeakProminence',1500000000);  
        % re-smooth data
        % Magn_each_blade = damping_re_smooth(Magn_each_blade,locs);
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
            % calculate damping ratio%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sol = roots([8,0,2,(f1 - f2) / freq_each_blade(locs(j))]); 
            damping_ratio(j) = sol(imag(sol) == 0 & sol >= 0 & sol <= 1);
            % plot peaks
            plot([freq_each_blade(locs(j)), freq_each_blade(locs(j))], [min(Magn_each_blade), max(Magn_each_blade)], 'r-');
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate damping ratios(single_degree_of_freedom_approximation)
if(method == "single_degree_of_freedom_approximation")    
    figure('units','normalized','outerposition',[0 0 1 1]);
    if (plot_struct == "all")    
        n_cols = 2;
        n_rows = ceil(n_blades / n_cols);
        num_blades = n_blades;
    end
    if (plot_struct == "one") 
         n_cols = 1;
         n_rows = 1;
         num_blades = 1;
    end    
    if (plot_struct == "four") 
         n_cols = 2;
         n_rows = ceil(4 / n_cols);
         num_blades = 4;
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
        % calculate damping ratio for each blade    
        Magn_each_blade = magn_i;
        freq_each_blade = freq_i;
        %find peaks(but not using smooth data)
        Magn_each_blade_smooth = smoothdata(Magn_each_blade,'loess',1000);           
        [pks,locs] = findpeaks(Magn_each_blade_smooth,'MinPeakProminence',10);  

        % deal each peak    
        damping_ratio = zeros(size(pks));
        for j = 1:length(pks)
            
    


            
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
% include using function damping_re_smooth(Magn_each_blade,locs)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




