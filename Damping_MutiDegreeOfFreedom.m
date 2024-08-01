function Damping_MutiDegreeOfFreedom(blade, EO,SetTag)
    fprintf('[****damping:muti degree of freedom approximation starts.****]\n');
    
    %% init params
    n_blades = length(blade);
    quality_factor_vec = [];

    %% get peaks_idx    
    peaks_idx_all = FindPeakAutomatic(blade);

    %% deal every blade
    for blade_idx = 1:n_blades
        %% init
        fprintf('blade:%d\n', blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 

        %% Create figure and set it up
        fig = figure('units', 'normalized', 'outerposition', [0 0.25 1 0.5]);
        sgtitle(sprintf('EO%d %s blade%d', EO,SetTag,blade_idx));
        set(gcf, 'WindowStyle', 'docked');
       
        %% reduce noise and downsample
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, 'movmean', 200);
        phase = smoothdata(phase,'movmean',200);

        %% subplot 1: Magnitude
        % subplot(2, 1, 1);
        plot(freq, magn, 'o', 'Color', [0.8, 0.9, 1.0], 'DisplayName', ...
            'Amplitude');  
        hold on;
        plot(freq, err, 'o', 'Color', [1.0, 0.8, 0.8], 'DisplayName', 'Error');  
        xlabel('Frequency(Hz)');
        ylabel('Normalized Amplitude');
        
        
        %% get peaks_idx              
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx);
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: can not find peaks for this blade, will skip');
            continue;
        end

        %% set weights
        weights_idx = Weight(magn, peaks_idx);

        %% set boundary
        boundary_idx = [weights_idx(1,1),weights_idx(end,end)];
        [freq_cut, magn_cut, weights_idx_cut, peaks_idx_cut] = SetBoundary(freq, ...
            magn, peaks_idx, weights_idx, boundary_idx);
        line([freq(boundary_idx(1)), freq(boundary_idx(1))], [0, max(magn)], ...
            'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Left Boundary');
        line([freq(boundary_idx(2)), freq(boundary_idx(2))], [0, max(magn)], ...
            'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Right Boundary');       

        %% least squared method
        params_m = LMAlgorithmMDOF(freq_cut, magn_cut, peaks_idx_cut, ...
            weights_idx_cut);
        plot(freq_cut, abs(ModelMDOF(params_m, freq_cut)), '--', 'Color', ...
            [0, 0.7, 0], 'DisplayName', 'Fitted Model');
        for i = 1:size(weights_idx_cut, 1)
            startIdx = weights_idx_cut(i, 1);
            endIdx = weights_idx_cut(i, 2); 
            plot(freq_cut(startIdx:endIdx), abs(ModelMDOF(params_m, ...
                freq_cut(startIdx:endIdx))),'-', 'Color', [0, 0.9, 0], ...
                'LineWidth', 6, 'DisplayName', 'Weight Area');
        end
        xlim([13680,13880]);                                             
        % ylim([]);
        hold off;       

        %% quality factor
        quality_factor = 0;
        quality_factor = abs(abs(ModelMDOF(params_m, freq_cut))-magn_cut);
        quality_factor = sum(quality_factor) /length(magn_cut);
        quality_factor = quality_factor / max(magn_cut);
        quality_factor_vec = [quality_factor_vec,quality_factor];
        fprintf("quality factor: %d",quality_factor);

        %% subplot 2: Phase
        % subplot(2, 1, 2); 
        % xlabel('Frequency(Hz)');
        % ylabel('Phase(Rad)');
        % hold on;
        % xlim([13650 13850]);                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 只适用于EO24!!
        % % from P_Phase
        % phase_from_data = rad2deg(phase);
        % plot(freq, phase, 'o', 'Color', [0.8, 0.9, 1.0], 'DisplayName', 'Phase from Data');

        %% save to file
        params_fitted{blade_idx} = reshape(params_m, 4, length(peaks_idx_cut)).';
        result_filename_excel = fullfile('Result', sprintf('EO%d_MDOF_%s.xlsx', ...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            EO,SetTag));
        T = table(params_fitted{blade_idx}(:, 1), params_fitted{blade_idx}(:, 2), ...
            params_fitted{blade_idx}(:, 3), params_fitted{blade_idx}(:, 4), ...
            'VariableNames', {'Frequency', 'D', 'r_re', 'r_im'});
        if isfile(result_filename_excel)
            existing_data = readtable(result_filename_excel);
            start_row = size(existing_data, 1) + 3;
        else
            start_row = 1;
        end
        writetable(T, result_filename_excel, 'Sheet', 1, 'Range', ...
            ['A' num2str(start_row)]);       

        % save graph
        graph_filename = fullfile('Graph', ...
            sprintf('EO%d_MDOF_blade%d_%s.png', EO, blade_idx,SetTag));              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        saveas(fig, graph_filename);    
    end
    % close all
    fprintf('[****damping:muti degree of freedom approximation finished.****]\n');
end
