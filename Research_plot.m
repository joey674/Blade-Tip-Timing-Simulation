function Research_plot(blades, EO)
    %% init params
    n_blades = length(blades);

    %% get peaks_idx    
    peaks_idx_all = FindPeakAutomatic(blades);

    %% Create figure and set it up
    fig = figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 600, 2000 * n_blades]); % Adjust the figure size
    sgtitle(sprintf('EO%d', EO));
    % set(gcf, 'WindowStyle', 'docked');

    %% deal with every blade
    for blade_idx = 1:n_blades
        subplot(n_blades,1,blade_idx); % Create a subplot for each blade

        %% init
        fprintf('blade:%d\n', blade_idx); 
        blade_data = blades{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 

        %% reduce noise and downsample
        magn = ReduceNoise(magn, freq, err);
        err = smoothdata(err, 'movmean', 200);

        %% plot the data
        plot(freq, magn);
         % xlabel('Frequency (Hz)');
        % ylabel('Magnitude');
        % title(sprintf('Blade %d', blade_idx));

        %% get peaks_idx              
        peaks_idx = [peaks_idx_all(blade_idx).peaks.idx];
        peaks_idx = sort(peaks_idx);
        if isempty(peaks_idx) 
            disp('ERROR: can not find peaks for this blade, will skip this blade.');
            continue;
        end

        %% mark the peaks
        hold on;
        plot(freq(peaks_idx), magn(peaks_idx), 'ro'); % Plot peaks in red 'o'
        hold off;
    end

    %% Save the figure
    saveas(fig, sprintf('peak_EO%d.png',EO));

    close all
end

