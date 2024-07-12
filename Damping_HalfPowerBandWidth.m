function Damping_HalfPowerBandWidth(blade,EO)
    fprintf('[**********damping:half power bandwidth method starts.**********]\n');
    
    %% init params
    n_blades = length(blade);

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

        %% Create figure and set it up
        fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
        sgtitle(sprintf('EO%d, blade%d', EO, blade_idx));
        set(gcf, 'WindowStyle', 'docked');

        
        %% reduce noise and downsample
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, 'movmean', 200);

        %% subplot: Magnitude
        plot(freq, magn, 'o', 'Color', [0.8, 0.9, 1.0], 'DisplayName', 'Magnitude');  
        hold on;
        plot(freq, err, 'o', 'Color', [1.0, 0.8, 0.8], 'DisplayName', 'Error');  
        xlabel('Frequency');
        ylabel('Magnitude');

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
            plot(freq(peaks_idx(j)), peaks_y(j), 'ro', 'DisplayName', 'Peak');
            plot([f1 f2], [peaks_y(j)/sqrt(2) peaks_y(j)/sqrt(2)], 'go', 'DisplayName', 'Half Power Point');
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
end
