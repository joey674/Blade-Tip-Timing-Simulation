function Damping_MutiDegreeOfFreedom(dataset_path,blade,EO)
    fprintf('[**********damping:muti degree of freedom approximation starts.**********]\n');
    
    %% init params
    n_blades = length(blade);
    EO_str = num2str(EO);
    filename = strrep(dataset_path, '.mat', ['_EO' EO_str '_PeaksIdx.mat']);
    if exist(filename,"file") == 2
       data = load(filename, 'peaks_idx_magn');
       % data = load(filename, 'peaks_idx_magn','weights_idx_magn');
    end 

    %% deal every blades
    for blade_idx = 1:n_blades 
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn]; 
        err  = [blade_data.err];          
        magn = Damping_NoiseFilter(magn);
        err = smoothdata(err,'movmean',100);

        %% get peaks_idx(load existed file or manually input)
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');
        title(sprintf('EO%d, blade%d', EO, blade_idx));xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
        hold on; 
        plot(freq,magn, 'Color', [0.6, 0.8, 1.0], 'DisplayName', 'Magnitude');plot(freq,err, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Error');     
        % legend;   
        peaks_idx = [];
        weights_idx = [];
        if exist(filename,"file") == 2% load existed preset file
            peaks_idx = data.peaks_idx_magn{blade_idx};
            % weights_idx = data.weights_idx_magn{blade_idx};
            for i = 1:length(peaks_idx)
                plot(freq(peaks_idx(i)), magn(peaks_idx(i)), 'bo', 'DisplayName', ['Peak idx ' num2str(peaks_idx(i))]);
                % plot(freq(weights_idx(i,1)), magn(weights_idx(i,1)), 'go', 'DisplayName', ['Left weight idx ' num2str(peaks_idx(i))]);
                % plot(freq(weights_idx(i,2)), magn(weights_idx(i,2)), 'go', 'DisplayName', ['Right weight idx ' num2str(peaks_idx(i))]);
            end
        end         
        if exist(filename,"file") ~= 2% manually input 
            while true
                [freq_get, ~] = ginput(1); % Get one point
                if isempty(freq_get)  % press enter to stop
                    peaks_idx_magn{blade_idx} = peaks_idx;
                    weights_idx_magn{blade_idx} = weights_idx;
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

                % title('Select left boundary');
                % [freq_get, ~] = ginput(1);
                % [~, left_weight_idx] = min(abs(freq - freq_get));
                % left_weight_idx = round(left_weight_idx);
                % plot(freq(left_weight_idx), magn(left_weight_idx), 'go', 'DisplayName', ['Left weight idx ' num2str(idx)]);
                % title('Select right boundary');
                % [freq_get, ~] = ginput(1);
                % [~, right_weight_idx] = min(abs(freq - freq_get));
                % right_weight_idx = round(right_weight_idx);
                % plot(freq(right_weight_idx), magn(right_weight_idx), 'go', 'DisplayName', ['Right weight idx ' num2str(idx)]);   
                % weights_idx = [weights_idx; [left_weight_idx, right_weight_idx]];
                % title(''); % Clear the title
            end
        end
        fprintf("peaks_idx:%d ",peaks_idx);fprintf("\n");    
        peaks_idx = sort(peaks_idx);

        %% set boundary and weights
        weights_idx = MDOF_Weight(magn,peaks_idx);
        boundary_idx = MDOF_Bound(magn,peaks_idx);
        [freq_cut,magn_cut,weights_idx_cut,peaks_idx_cut] = MDOF_SetBoundary(freq,magn,peaks_idx,weights_idx,boundary_idx);
        line([freq(boundary_idx(1)),freq(boundary_idx(1))],[0,max(magn)],'Color','red','LineStyle','--');
        line([freq(boundary_idx(2)),freq(boundary_idx(2))],[0,max(magn)],'Color','red','LineStyle','--');       

        %% least squared method
        params_m = MDOF_LMAlgorithm(freq_cut,magn_cut,peaks_idx_cut,weights_idx_cut);
        plot(freq_cut, MDOF_Model(params_m,freq_cut),'--','Color',[0,0.9,0], 'DisplayName', 'Fitted Model');
        for i = 1:size(weights_idx_cut, 1)
            startIdx = weights_idx_cut(i, 1);
            endIdx = weights_idx_cut(i, 2);
            plot(freq_cut(startIdx:endIdx), MDOF_Model(params_m,freq_cut(startIdx:endIdx)),'Color',[0,0.7,0]);
        end
        hold off;       

        %% save to file
        params_fitted{blade_idx} = reshape(params_m,4, length(peaks_idx_cut)).';
        excitate_freq {blade_idx} = params_fitted{blade_idx}(:,1).';
        damping_ratios{blade_idx} = params_fitted{blade_idx}(:,2).';
        processed_magn{blade_idx} = magn;% for DL  
        
        %% save graph
        graphname = sprintf('graph\\EO%d_blade%d_LM_NonLinWeight.png', EO, blade_idx);% save graph
        saveas(gcf, graphname);        
    end
    if exist(filename,"file") ~= 2% save to file
        save(filename, 'peaks_idx_magn','weights_idx_magn','processed_magn','damping_ratios','excitate_freq');
        save(filename, 'peaks_idx_magn','processed_magn','damping_ratios','excitate_freq');
    end


    fprintf('[**********damping:muti degree of freedom approximation finished.**********]\n');
end

