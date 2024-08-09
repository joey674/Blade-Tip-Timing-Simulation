function Damping_MutiDegreeOfFreedom(blade, EO, Tag)
    fprintf('[****damping:multi degree of freedom approximation starts.****]\n');
    
    %% Initialize parameters
    %{
        x_left/x_right: used for the axis of the plot;
        normalized_factor: used to plot;
            stands for the max of all the Amplitude; 
            or can use the max of the specific blade below: normalized_factor = max(magn)+10;
    %} 
    quality_factor_vec = [];
    if EO == 24
        x_left = 13680;
        x_right = 13880;
        normalized_factor = 120;
    elseif EO == 20
         x_left = 9420;
        x_right = 9730;
        normalized_factor = 55;
    elseif EO == 8
         x_left = 3770;
        x_right = 3890;
        normalized_factor = 65;
    end

    %% interpolation to the same distance
    %{
        interpolate the data to the same freq density
    %}
    blade = Interpolate(blade);

    %% Get peaks indices
    %{
        peaks_idx_all: contains peaks idx for all blades
    %}
    peaks_idx_all = FindPeakAutomatic(blade);

    %% Process each blade
    for blade_idx = 1:12
        %% Initialize
        fprintf('blade:%d\n', blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 
       
        %% Reduce noise 
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
        
        %% Get peaks indices
        %{
            press enter to the next blade plot;
            "[fig, peaks_idx] = FindPeakManual(fig, freq, magn,
            peaks_idx)" can be muted seperatly to make a fast browse 
        %}
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx);
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: cannot find peaks for this blade, will skip');
            continue;
        end

        %% Set weights
        weights_idx = Weight(magn, peaks_idx);

        %% Set boundary
        boundary_idx = [weights_idx(1, 1), weights_idx(end, end)];
        [freq_cut, magn_cut, weights_idx_cut, peaks_idx_cut] = SetBoundary(freq, magn, peaks_idx, weights_idx, boundary_idx);
        line([freq(boundary_idx(1)), freq(boundary_idx(1))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Left Boundary');
        line([freq(boundary_idx(2)), freq(boundary_idx(2))], [0, max(magn)], 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'Right Boundary');       

        %% Least squares method
        params_m = LMAlgorithmMDOF(freq_cut, magn_cut, peaks_idx_cut, weights_idx_cut);
        plot(freq_cut, abs(ModelMDOF(params_m, freq_cut)), '--', 'Color', [0, 0.7, 0], 'DisplayName', 'Fitted Model');
        for i = 1:size(weights_idx_cut, 1)
            startIdx = weights_idx_cut(i, 1);
            endIdx = weights_idx_cut(i, 2); 
            plot(freq_cut(startIdx:endIdx), abs(ModelMDOF(params_m, freq_cut(startIdx:endIdx))), '-', 'Color', [0, 0.9, 0], 'LineWidth', 6, 'DisplayName', 'Weight Area');
        end
        hold off;       

        %% Calculate Quality factor
        quality_factor = abs(abs(ModelMDOF(params_m, freq_cut)) - magn_cut);
        quality_factor = sum(quality_factor) / length(magn_cut);
        quality_factor = 1-quality_factor / max(magn_cut);
        quality_factor_vec = [quality_factor_vec, quality_factor];
        fprintf("quality factor: %.2f%%\n", quality_factor * 100);
        % add to the graph
        x_limits = xlim;
        y_limits = ylim;
        x_pos = x_limits(2) - 0.05 * diff(x_limits); 
        y_pos = y_limits(2) - 0.05 * diff(y_limits); 
        text(x_pos, y_pos, sprintf('Quality Factor: %.2f%%', quality_factor * 100), ...
             'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');

        %% Save results to file
        params_fitted{blade_idx} = reshape(params_m, 4, length(peaks_idx_cut)).';
        result_filename_excel = fullfile('Result', sprintf('EO%d_MDOF_%s.xlsx', EO, Tag));
        T = table(params_fitted{blade_idx}(:, 1), params_fitted{blade_idx}(:, 2), params_fitted{blade_idx}(:, 3), params_fitted{blade_idx}(:, 4), ...
            'VariableNames', {'Frequency', 'D', 'r_re', 'r_im'});
        if isfile(result_filename_excel)
            existing_data = readtable(result_filename_excel);
            start_row = size(existing_data, 1) + 3;
        else
            start_row = 1;
        end
        writetable(T, result_filename_excel, 'Sheet', 1, 'Range', ['A' num2str(start_row)]);       

        %% Save graph
        graph_filename = fullfile('Graph', sprintf('EO%d_MDOF_blade%d_%s.png', EO, blade_idx, Tag));
        saveas(fig, graph_filename);    
    end
    close all
    fprintf('[****damping:multi degree of freedom approximation finished.****]\n');
end

