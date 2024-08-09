function Damping_SingleDegreeOfFreedom(blade,EO,Tag)
    fprintf('[**********damping:SDOF method starts.**********]\n');
    
    %% Initialize parameters
    n_blades = length(blade);
    quality_factor_vec = [];
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
    %{
        interpolate the data to the same freq density
    %}
    blade = Interpolate(blade);

    %% get peaks_idx    
    peaks_idx_all = FindPeakAutomatic(blade);

    %% deal every blades
    for blade_idx = 1:12
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 

        %% Reduce noise and downsample
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, 'movmean', 200);
        phase = smoothdata(phase, 'movmean', 200);
        
        %% Normalize
        normalized_factor = max(magn)+10;
        magn = magn / normalized_factor;
        err = err / normalized_factor;

        %% Create figure and set it up
        fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 400]);
        sgtitle(sprintf('EO%d %s blade%d', EO, Tag, blade_idx));
        set(fig, 'WindowStyle', 'normal');
        plot(freq, magn, 'o', 'Color', [0.8, 0.9, 1.0], 'DisplayName', 'Amplitude');  
        hold on;
        plot(freq, err, 'o', 'Color', [1.0, 0.8, 0.8], 'DisplayName', 'Error');  
        xlim([x_left,x_right]);
        ylim([0, 1]);
        xlabel('Frequency (Hz)');
        ylabel('Normalized Amplitude');

        %% get peaks_idx              
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx);
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: can not find peaks for this blade, will skip this blade.');
            continue;
        end        
        
        %% fit for every mode
        damping_ratio = zeros(size(peaks_idx));
        excitate_freq = zeros(size(peaks_idx));
        for j = 1:length(peaks_idx)
            %% set weights
            weight_idx = Weight(magn, peaks_idx(j));   
    
            %% set boundary
            boundary_idx = [weight_idx(1),weight_idx(2)];
            [freq_cut, magn_cut, weights_idx_cut, peaks_idx_cut] = SetBoundary(freq, magn, peaks_idx(j), weight_idx, boundary_idx);
            line([freq(boundary_idx(1)), freq(boundary_idx(1))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--');
            line([freq(boundary_idx(2)), freq(boundary_idx(2))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--');       
    
            %% least squared method
            params_m = LMAlgorithmSDOF(freq_cut, magn_cut, peaks_idx_cut, weights_idx_cut);
            plot(freq_cut, abs(ModelSDOF(params_m, freq_cut)), '--', 'Color', [0, 0.7, 0]);
            for i = 1:size(weights_idx_cut, 1)
                startIdx = weights_idx_cut(i, 1);
                endIdx = weights_idx_cut(i, 2); 
                plot(freq_cut(startIdx:endIdx), abs(ModelSDOF(params_m, freq_cut(startIdx:endIdx))), '-', 'Color', [0, 0.9, 0], 'LineWidth', 6, 'DisplayName', 'Weight Area');
            end
            excitate_freq(j) = params_m(1);
            damping_ratio(j) = params_m(2);

            %% Quality factor
            quality_factor = abs(abs(ModelSDOF(params_m, freq_cut)) - magn_cut);
            quality_factor = sum(quality_factor) / length(magn_cut);
            quality_factor = 1-quality_factor / max(magn_cut);
            quality_factor_vec = [quality_factor_vec, quality_factor];
            fprintf("quality factor: %.2f%%\n", quality_factor * 100);
          
        end
        x_limits = xlim;
        y_limits = ylim;

        x_pos = x_limits(2) - 0.05 * diff(x_limits);
        y_pos = y_limits(2) - 0.05 * diff(y_limits);
        
        quality_factors_str = sprintf('Quality Factors: %s', ...
        strjoin(arrayfun(@(i, qf) sprintf('QF %d: %.2f%%', i, qf * 100), ...
        1:length(quality_factor_vec), quality_factor_vec, 'UniformOutput', false), '; '));
        
        % Add the concatenated text to the plot
        text(x_pos, y_pos, quality_factors_str, ...
            'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');
        quality_factor_vec = [];
      
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
    close all
    fprintf('[**********damping:SDOF Method finished.**********]\n');
end
