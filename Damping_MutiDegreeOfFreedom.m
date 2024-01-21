function Damping_MutiDegreeOfFreedom(dataset_path,blade,EO)
    fprintf('[**********damping:muti degree of freedom approximation starts.**********]\n');
    
    %% init params
    n_blades = length(blade);
    EO_str = num2str(EO);
    filename = strrep(dataset_path, '.mat', ['_EO' EO_str '_PeaksIdx.mat']);
    if exist(filename,"file") == 2
       data = load(filename, 'peaks_idx_magn');
    end 

    %% deal every blades
    for blade_idx = 1:n_blades 
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};% import data;
        freq = [blade_data.freq];
        magn = [blade_data.magn]; 
        err  = [blade_data.err];          
        magn = Damping_NoiseFilter(magn);% smooth data 
        err = smoothdata(err,'movmean',100);

        %% get peaks_idx(load existed file or manually input)
        figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title("blade",blade_idx);xlabel('Frequency (Hz)');ylabel('Magnitude (mm)');
        hold on; 
        plot(freq,magn, 'Color', [0.6, 0.8, 1.0], 'DisplayName', 'Magnitude');plot(freq,err, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Error');     
        legend;   
        peaks_idx = [];
        if exist(filename,"file") == 2% load existed preset file
            peaks_idx = data.peaks_idx_magn{blade_idx};
            for i = 1:length(peaks_idx)
                plot(freq(peaks_idx(i)), magn(peaks_idx(i)), 'bo', 'DisplayName', ['Peak idx ' num2str(peaks_idx(i))]);
            end
        end         
        if exist(filename,"file") ~= 2% manually input 
            while true
                [freq_get, ~] = ginput(1); % Get one point
                if isempty(freq_get)  % press enter to stop
                    peaks_idx_magn{blade_idx} = peaks_idx;
                    break;
                end
                [~, idx] = min(abs(freq - freq_get));
                idx = round(idx);
                searchRange = 1000;
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

        %% set boundary and weights
        weights_idx = MDOF_Weight(magn,peaks_idx);
        boundary_idx = MDOF_Bound(magn,peaks_idx);
        [freq_cut,magn_cut,weights_idx_cut,peaks_idx_cut] = MDOF_SetBoundary(freq,magn,peaks_idx,weights_idx,boundary_idx);
        % plot 
        line([freq(boundary_idx(1)),freq(boundary_idx(1))],[0,max(magn)],'Color','red','LineStyle','--','DisplayName', 'left boundary');
        line([freq(boundary_idx(2)),freq(boundary_idx(2))],[0,max(magn)],'Color','red','LineStyle','--','DisplayName', 'right boundary');       
        % plot(freq_cut(weights_idx_cut(:,1)+100), magn_cut(weights_idx_cut(:,1)), 'ro','DisplayName', ['left weight idx']);
        % plot(freq_cut(weights_idx_cut(:,2)-100), magn_cut(weights_idx_cut(:,2)), 'yo','DisplayName', ['left weight idx']);
        
        %% least squared method
        params_m = MDOF_LMAlgorithm(freq_cut,magn_cut,peaks_idx_cut,weights_idx_cut);
        plot(freq_cut, MDOF_Model(params_m,freq_cut),'--','Color',[0,0.9,0], 'DisplayName', 'Fitted Model');
        for i = 1:size(weights_idx_cut, 1)
            startIdx = weights_idx_cut(i, 1);
            endIdx = weights_idx_cut(i, 2);
            plot(freq_cut(startIdx:endIdx), MDOF_Model(params_m,freq_cut(startIdx:endIdx)),'Color',[0,0.7,0],'DisplayName', ['Weighted Range ' num2str(i)]);
        end
        hold off;       

        %% save to file
        params_fitted{blade_idx} = reshape(params_m,4, length(peaks_idx_cut)).';
        excitate_freq {blade_idx} = params_fitted{blade_idx}(:,1).';
        damping_ratios{blade_idx} = params_fitted{blade_idx}(:,2).';
        processed_magn{blade_idx} = magn;% for DL  
        
        %% save graph
        % graphname = sprintf('graph\\EO%d_blade%d.png', EO, blade_idx);% save graph
        % saveas(gcf, graphname);        
    end
    if exist(filename,"file") ~= 2% save to file
        save(filename, 'peaks_idx_magn','processed_magn','damping_ratios','excitate_freq');
    end


    fprintf('[**********damping:muti degree of freedom approximation finished.**********]\n');
end

