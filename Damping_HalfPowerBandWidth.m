function Damping_HalfPowerBandWidth(blade,EO,Tag)
    fprintf('[**********damping:half power bandwidth method starts.**********]\n');
      
    %% Initialize parameters
    n_blades = length(blade);
    if EO == 24
        x_left = 13680;
        x_right = 13880;
    elseif EO == 20
         x_left = 9420;
        x_right = 9730;
    elseif EO == 8
         x_left = 3770;
        x_right = 3890;
    end

    %% interpolation to the same distance
    blade = Interpolate(blade);

    %% get peaks_idx    
    peaks_idx_all = FindPeakAutomatic(blade);

    %% deal every blades
    for blade_idx = 1:n_blades 
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 

        
        %% reduce noise and downsample
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, 'movmean', 200);

        %% Normalize
        normalized_factor = max(magn)+10;
        magn = magn / normalized_factor;
        err = err / normalized_factor;


        %% Create figure and set it up
        fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 400]);
        sgtitle(sprintf('EO%d %s blade%d', EO,Tag,blade_idx));
        set(fig, 'WindowStyle', 'normal');
        plot(freq, magn, 'o', 'Color', [0.8, 0.9, 1.0], 'DisplayName', ...
            'Amplitude');  
        hold on;
        plot(freq, err, 'o', 'Color', [1.0, 0.8, 0.8], 'DisplayName', 'Error');  
        xlim([x_left,x_right]);
        ylim([0,1]);
        xlabel('Frequency(Hz)');
        ylabel('Normalized Amplitude');

        %% get peaks_idx              
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx);
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: can not find peaks for this blade, will skip this blade.');
            continue;
        end        
        
        %% calculate damping ratio
        damping_ratio = zeros(size(peaks_idx));
        peaks_y = magn(peaks_idx);
        for j = 1:length(peaks_idx)
            % find half bandwith point 
            f1 = freq(find(magn(1:peaks_idx(j)) <= peaks_y(j)/sqrt(2), 1, 'last'));
            f2 = freq(find(magn(peaks_idx(j):end) <= peaks_y(j)/sqrt(2), 1, 'first') + peaks_idx(j) - 1);
            % check if f1, f2 are found
            if isempty(f1)
                f1 = freq(1);
            end
            if isempty(f2)
                f2 = freq(end);
            end
            % calculate damping ratio
            damping_ratio(j) = (f2 - f1) / (2 * freq(peaks_idx(j)));
            
            % plot peak and half power points
            plot(freq(peaks_idx(j)), peaks_y(j), 'o', 'Color', 'r', 'MarkerFaceColor', 'r');
            plot([f1 f2], [peaks_y(j)/sqrt(2) peaks_y(j)/sqrt(2)], 'o', 'Color', 'g', 'MarkerFaceColor', 'g');
        end    

        %% save to file 
        damping_ratios{blade_idx} = damping_ratio';
        excitate_freq{blade_idx} = freq(peaks_idx)';
        result_filename_excel = fullfile('Result', sprintf('EO%d_HPBW.xlsx', EO));
        T = table(excitate_freq{blade_idx},damping_ratios{blade_idx}, ...
            'VariableNames', {'Frequency', 'D'});
        if isfile(result_filename_excel)
            existing_data = readtable(result_filename_excel);
            start_row = size(existing_data, 1) + 3;
        else
            start_row = 1;
        end
        writetable(T, result_filename_excel, 'Sheet', 1, 'Range', ['A' num2str(start_row)]);       
        disp(['Data saved to ', result_filename_excel]);

        %% save graph
        graph_filename = fullfile('Graph', sprintf('EO%d_HPBW_blade%d.png', EO, blade_idx));
        saveas(fig, graph_filename);    
    end        
    fprintf('[**********damping:half power bandwidth finished.**********]\n');
    close all;
end
