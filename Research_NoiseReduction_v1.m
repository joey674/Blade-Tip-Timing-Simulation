function Research_NoiseReduction_v1(blade,EO)
    %% init params
    n_blades = length(blade);
    peakidx_filename = fullfile('PeaksIdx', sprintf('EO%d_PeaksIdx.mat', EO));
    result_filename =  fullfile('Result', sprintf('EO%d_MDOF.mat', EO));
    if exist(peakidx_filename,"file") == 2
       data = load(peakidx_filename, 'peaks_idx_magn');
    end 

        %% deal every blades
    for blade_idx = 1:1
        %% init
        fprintf('blade:%d\n',blade_idx); 
        blade_data = blade{blade_idx};
    
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err  = [blade_data.err]; 
        % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title('origin');
        % plot(magn, 'o');hold on;
        % % plot(err,'o');
        % hold off;
    
        %% 查找密度范围
        density_levels =        [1,   0.4, 0.1, 0.05, 0.01]; 
        downsample_factors =    [10 , 15,  20,  40,   50  ]; 
        [err_density, err_values] = ksdensity(err);

        % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title('error density');
        % plot(err_values,err_density, 'o', 'DisplayName', 'error density');hold on;
        % hold off;

        [max_density, max_index] = max(err_density);
        err_value_lower_last = err_values(max_index);
        err_value_upper_last = err_values(max_index);
        downsampled_freq = freq;
        downsampled_magn = magn;      
        for i = 1:length(density_levels)
                target_density = density_levels(i) * max_density;
                left_index = max_index;
                while left_index > 1 && err_density(left_index) >= target_density
                    left_index = left_index - 1;
                end
                right_index = max_index;
                while right_index < length(err_density) && err_density(right_index) >= target_density
                    right_index = right_index + 1;
                end
                err_value_lower_current = err_values(left_index);
                err_value_upper_current = err_values(right_index);

                %% 
                indices_downsampled = [ find(downsampled_magn > err_value_lower_current & downsampled_magn < err_value_lower_last),...
                                        find(downsampled_magn < err_value_upper_current & downsampled_magn > err_value_upper_last)];
                indices_downsampled = custom_downsample(indices_downsampled, downsample_factors(i));
                indices_outside = [find(downsampled_magn < err_value_lower_current),...
                                   find(downsampled_magn > err_value_upper_current),...
                                   find(downsampled_magn > err_value_lower_last & downsampled_magn < err_value_upper_last)];

                indices_all = [indices_downsampled,indices_outside];
                
                downsampled_magn = downsampled_magn(sort(indices_all)); 
                downsampled_freq = downsampled_freq(sort(indices_all)); 
                
                err_value_lower_last = err_value_lower_current;
                err_value_upper_last = err_value_upper_current;
                % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title('freq response after downsample');
                % plot(downsampled_freq,downsampled_magn, 'o', 'DisplayName', 'Downsampled Magnitude Data');hold on;
                % hold off;
        end
        % figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title('freq response after downsample');
        % % plot(freq,magn, 'o');hold on;
        % plot(freq,smoothdata(magn,'movmean',100), 'LineWidth',2, 'DisplayName', 'Movmean');hold on;
        % plot(downsampled_freq,smoothdata(downsampled_magn,'movmean',100),'LineWidth',2,'DisplayName', 'movmean with downsample');hold on;
        % legend;
        % hold off;

        figure('units', 'normalized', 'outerposition', [0 0 1 1]);set(gcf, 'WindowStyle', 'docked');title('method comparsion');
        % plot(freq,smoothdata(magn,'movmean',100), 'LineWidth',2, 'DisplayName', 'Movmean');hold on;
        % plot(downsampled_freq,smoothdata(downsampled_magn,'movmean',200),'LineWidth',1,'DisplayName', 'movmean');hold on;
        % plot(downsampled_freq,smoothdata(downsampled_magn,'rlowess',200),'LineWidth',1,'DisplayName', 'rlowess');hold on;
        % plot(downsampled_freq,smoothdata(downsampled_magn,'sgolay',200),'LineWidth',1,'DisplayName', 'sgoly');hold on;
        plot(downsampled_freq,smoothdata(downsampled_magn,'gaussian',200),'LineWidth',1,'DisplayName', 'gaussian');hold on;
        legend;
        hold off;
    end
end