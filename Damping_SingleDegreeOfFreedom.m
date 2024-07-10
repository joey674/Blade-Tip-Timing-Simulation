function Damping_SingleDegreeOfFreedom(blade,EO)
    fprintf('[**********damping:SDOF method starts.**********]\n');
    
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
        legend;

        %% get peaks_idx              
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        % [fig, peaks_idx] = MDOF_FindPeakManual(fig, freq, magn, peaks_idx);
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: can not find peaks for this blade, will skip this blade.');
            continue;
        end        
        
        damping_ratio = zeros(size(peaks_idx));
        excitate_freq = zeros(size(peaks_idx));
        for j = 1:length(peaks_idx)
            %% set weights
            weight_idx = Weight(magn, peaks_idx(j));   
    
            %% set boundary
            boundary_idx = [weight_idx(1),weight_idx(2)];
            [freq_cut, magn_cut, weights_idx_cut, peaks_idx_cut] = SetBoundary(freq, magn, peaks_idx(j), weight_idx, boundary_idx);
            line([freq(boundary_idx(1)), freq(boundary_idx(1))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Left Boundary');
            line([freq(boundary_idx(2)), freq(boundary_idx(2))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Right Boundary');       
    
            %% least squared method
            params_m = LMAlgorithmSDOF(freq_cut, magn_cut, peaks_idx_cut, weights_idx_cut);
            plot(freq_cut, abs(ModelSDOF(params_m, freq_cut)), '--', 'Color', [0, 0.7, 0], 'DisplayName', 'Fitted Model');
            for i = 1:size(weights_idx_cut, 1)
                startIdx = weights_idx_cut(i, 1);
                endIdx = weights_idx_cut(i, 2); 
                plot(freq_cut(startIdx:endIdx), abs(ModelSDOF(params_m, freq_cut(startIdx:endIdx))), 'o', 'Color', [0, 0.9, 0], 'DisplayName', 'Weight Area');
            end
            legend;    
            excitate_freq(j) = params_m(1);
            damping_ratio(j) = params_m(2);
        end

        %% save to file 
        damping_ratios{blade_idx} = damping_ratio';
        excitate_freqs{blade_idx} = freq(peaks_idx)';
        result_filename_excel = fullfile('Result', sprintf('EO%d_SDOF.xlsx', EO));
        T = table(excitate_freqs{blade_idx},damping_ratios{blade_idx}, ...
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
        graph_filename = fullfile('Graph', sprintf('EO%d_SDOF_blade%d.png', EO, blade_idx));
        saveas(fig, graph_filename);    
    end        
    fprintf('[**********damping:SDOF Method finished.**********]\n');
end
