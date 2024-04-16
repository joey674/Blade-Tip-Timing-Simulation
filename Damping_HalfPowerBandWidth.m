function Damping_HalfPowerBandWidth(blade,EO)
    fprintf('[**********damping:half power bandwidth method starts.**********]\n');
    
    %% init params
    n_blades = length(blade);
    peakidx_filename = fullfile('PeaksIdx', sprintf('EO%d_PeaksIdx.mat', EO));
    result_filename =  fullfile('Result', sprintf('EO%d_HPBW.mat', EO));
    if exist(peakidx_filename,"file") == 2
       data = load(peakidx_filename, 'peaks_idx_magn');
    end 

    %% deal every blades
    for blade_idx = 1:n_blades 
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 
        
        %% plot smoothed magn
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);subplot(2,1,1);set(gcf, 'WindowStyle', 'docked');
        title(sprintf('EO%d, blade%d', EO, blade_idx));xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
        hold on;
        legend;
        magn = smoothdata(magn,"movmean",100);
        plot(freq, magn,'Color', [0.7, 0.8, 1.0],'DisplayName', 'Magnitude after movmean');  

        %% get peaks_idx(load existed file)               
        peaks_idx = [];
        if exist(peakidx_filename,"file") == 2% load existed preset file
            peaks_idx = data.peaks_idx_magn{blade_idx};
        else 
            err("need to select peak in MDOF first");
        end

        %% calculate damping ratio
        damping_ratio = zeros(size(peaks_idx));
        peaks_y = magn(peaks_idx);
        for j = 1:length(peaks_idx)
            % find half bandwith point 
            f1 = freq(find(magn(1:peaks_idx(j)) <= peaks_y(j)/sqrt(2), 1, 'last'));
            f2 = freq(find(magn(peaks_idx(j):end) <= peaks_y(j)/sqrt(2), 1, 'first') + peaks_idx(j) - 1);
            %check if f1,f2 are found
            if isempty(f1)
                f1 = freq(1);
            end
            if isempty(f2)
                f2 = freq(end);
            end
            % calculate damping ratio
            damping_ratio(j) = (f2 - f1) / (2 * freq(peaks_idx(j)));
        end    

        %% save to file 
        damping_ratios{blade_idx} = damping_ratio;
        excitate_freq{blade_idx} = freq(peaks_idx);
    end        
    save(result_filename,'damping_ratios','excitate_freq');

    fprintf('[**********damping:half power bandwidth finished.**********]\n');
end

