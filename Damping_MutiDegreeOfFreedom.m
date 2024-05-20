function Damping_MutiDegreeOfFreedom(blade,EO)
    fprintf('[**********damping:muti degree of freedom approximation starts.**********]\n');
    
    %% init params
    n_blades = length(blade);
    peakidx_filename = fullfile('PeaksIdx', sprintf('EO%d_PeaksIdx.mat', EO));
    result_filename =  fullfile('Result', sprintf('EO%d_MDOF.mat', EO));
    if exist(peakidx_filename,"file") == 2
       data = load(peakidx_filename, 'peaks_idx_magn');
    end 

    %% deal every blades
    % for blade_idx = 1:n_blades
    for blade_idx = 1:1
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 
       
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);subplot(2,1,1);set(gcf, 'WindowStyle', 'docked');
        title(sprintf('EO%d, blade%d', EO, blade_idx));xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
        hold on;
        legend;
       
        %% noise reduce
        magn = Damping_NoiseFilter(magn);
        err = Damping_NoiseFilter(err);
        % phase = smoothdata(phase,'movmean',10);
        
        plot(freq, magn,'Color', [0.7, 0.8, 1.0],'DisplayName', 'Magnitude after movmean');  

        %% get peaks_idx(load existed file or manually input)               
        peaks_idx = [];
        if exist(peakidx_filename,"file") == 2% load existed preset file
            peaks_idx = data.peaks_idx_magn{blade_idx};
            % for i = 1:length(peaks_idx)
            %     plot(freq(peaks_idx(i)), magn(peaks_idx(i)), 'bo', 'DisplayName', ['Peak idx ' num2str(peaks_idx(i))]);
            % end
        end         
        if exist(peakidx_filename,"file") ~= 2% manually input 
            % error('manually input is unimplemented');
            while true
                [freq_get, ~] = ginput(1); % Get one point
                if isempty(freq_get)  % press enter to stop
                    break;
                end
                [~, idx] = min(abs(freq - freq_get));
                idx = round(idx);
                searchRange = 300;
                startIndex = max(1, idx - searchRange);
                endIndex = min(length(magn), idx + searchRange);
                [~, maxIndex] = max(magn(startIndex:endIndex));
                idx = startIndex + maxIndex - 1;
                plot(freq(idx), magn(idx), 'bo', 'DisplayName', ['Peak idx ' num2str(idx)]);
                peaks_idx = [peaks_idx,idx];            
            end
        end
        fprintf("peaks_idx:%d ",peaks_idx);fprintf("\n");    
        peaks_idx = sort(peaks_idx);

        %% set weights
        weights_idx = MDOF_Weight(magn,peaks_idx);

        %% set boundary
        boundary_idx = MDOF_Bound(magn,peaks_idx);
    
        %% cut unneeded part to plot and simulate
        [freq_cut,magn_cut,weights_idx_cut,peaks_idx_cut] = MDOF_SetBoundary(freq,magn,peaks_idx,weights_idx,boundary_idx);
        line([freq(boundary_idx(1)),freq(boundary_idx(1))],[0,max(magn)],'Color','red','LineStyle','--','DisplayName','Left Boundary');
        line([freq(boundary_idx(2)),freq(boundary_idx(2))],[0,max(magn)],'Color','red','LineStyle','--','DisplayName','Right Boundary');       

        %% least squared method
        params_m = MDOF_LMAlgorithm(freq_cut,magn_cut,peaks_idx_cut,weights_idx_cut);
        plot(freq_cut, abs(MDOF_Model(params_m,freq_cut)),'--','Color',[0,0.7,0], 'DisplayName', 'Fitted Model');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:size(weights_idx_cut, 1)
            startIdx = weights_idx_cut(i, 1);
            endIdx = weights_idx_cut(i, 2); 
            plot(freq_cut(startIdx:endIdx), abs(MDOF_Model(params_m,freq_cut(startIdx:endIdx))),'Color',[0,0.9,0],'DisplayName','Weight Area');
        end
        hold off;       

        %% plot for phase
        subplot(2, 1, 2); 
        title(sprintf('Phase'));
        xlabel('Frequency (Hz)');ylabel('Phase (degrees)');
        hold on;
        % from P_Phase
        phase_from_data = phase;
        phase_from_data = rad2deg(phase);
        % phase_from_data = unwrap(phase_from_data);
        plot(freq, phase_from_data, 'Color', [0.8, 0.9, 1.0]);
        plot(freq, smoothdata(phase_from_data,"movmean",100), 'Color', [0.7, 0.8, 1.0]);
        % from model
        phase_from_model = atan(imag(MDOF_Model(params_m, freq_cut)) ./ real(MDOF_Model(params_m, freq_cut)));
        phase_from_model = rad2deg(phase_from_model);
        % phase_from_model = unwrap(phase_from_model);
        plot(freq_cut, phase_from_model, 'Color', [0, 0.9, 0], 'DisplayName', 'Phase');
        line([freq(boundary_idx(1)),freq(boundary_idx(1))],[-150,150],'Color','red','LineStyle','--');
        line([freq(boundary_idx(2)),freq(boundary_idx(2))],[-150,150],'Color','red','LineStyle','--');   
        hold off;

        %% save to file
        params_fitted{blade_idx} = reshape(params_m,3, length(peaks_idx_cut)).';
        excitate_freq{blade_idx} = params_fitted{blade_idx}(:,1).';
        damping_ratios{blade_idx} = params_fitted{blade_idx}(:,2).'; 
        peaks_idx_magn{blade_idx} = peaks_idx;
        excitate_phase{blade_idx} = phase_from_model(peaks_idx_cut);
        
        %% save graph
        % graphname = sprintf('graph\\EO%d_blade%d_LM_NonLinWeight_WithPhase_ReDenoise.png', EO, blade_idx);% save graph
        % saveas(gcf, graphname);        
    end
    save(peakidx_filename, 'peaks_idx_magn');
    save(result_filename,'excitate_phase','damping_ratios','excitate_freq');

    fprintf('[**********damping:muti degree of freedom approximation finished.**********]\n');
end

