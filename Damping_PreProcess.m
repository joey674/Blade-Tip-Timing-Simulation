function blade = Damping_PreProcess(EO,data,n_blades)
    blade = cell(1,length(EO));
    for i = 1:length(EO)
        Freq = data.mean_RPM * EO(i) / 60; % 数据Magn/Err是一样的，只是Freq是乘以不同的EO     
        Magn = cell(1, n_blades); 
        Phase = cell(1, n_blades); 
        Err = cell(1, n_blades);         
        for j = 1:n_blades
            % EO 24 16 
            if length(EO) == 1
                Magn{j} = data.P_Magn{j,EO(1)};
                Phase{j} = data.P_Phase{j,EO(1)}; 
                Err{j} = data.Fit_Error{j,EO(1)};               
            end
            % EO 8/20
            if length(EO) == 2 
                Magn{j} = data.P_Magn{j,EO(2),EO(1)}(:,i);  
                Phase{j} = data.P_Phase{j,EO(2),EO(1)}(:,i);  
                Err{j} = data.Fit_Error{j,EO(2),EO(1)}; 
            end                       
        end 

        % pre process of data
        blade{i} = cell(1, n_blades);   
        % find max length
        max_length = max([length(Freq), cellfun(@length, Magn)]);   
        for k = 1:n_blades
            % assign Magn{i} and Freq to the structure
            blade{i}{k} = repmat(struct('freq', NaN, 'magn', NaN, 'err', NaN ), 1, max_length);       
            len_freq = length(Freq);
            len_magn = length(Magn{k});
            for j = 1:max_length
                if j <= len_freq
                    blade{i}{k}(j).freq = Freq(j);
                end
                if j <= len_magn
                    blade{i}{k}(j).magn = Magn{k}(j);
                    blade{i}{k}(j).phase = Phase{k}(j);
                    blade{i}{k}(j).err = Err{k}(j);
                end
            end    
            % sort the structure array by freq
            [~, idx_sort] = sort([blade{i}{k}.freq]);
            blade{i}{k} = blade{i}{k}(idx_sort);   
            % remove elements contain nan
            idx_remove = isnan([blade{i}{k}.freq]) | isnan([blade{i}{k}.magn]);                   
            blade{i}{k}(idx_remove) = [];   
        end
        
    end
end

