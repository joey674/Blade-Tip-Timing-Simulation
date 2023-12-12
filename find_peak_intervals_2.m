function peak_intervals = find_peak_intervals_2(n_blades, n_modes, blade)
    peak_interval_width = 10; 

    for i = 1:n_blades
        lens(i) = length([blade{i}.magn]);
    end
    magn_all = NaN(max(lens),n_blades);
    for i = 1:n_blades
        magn_all(1:length([blade{i}.magn]),i) = [blade{i}.magn];
    end
    magn_max=max(magn_all,[],2,"omitnan");

    % init
    peak_intervals = zeros(2, n_modes);
    % find peak
    [pks, locs] = findpeaks(magn_max, 'MinPeakProminence', 5); 
    % create intervals
    freq = [blade{1}.magn];
    for i = 1:min(n_modes, length(pks))
        lower_bound = max(freq(locs(i)) - peak_interval_width/2, min(freq));
        upper_bound = min(freq(locs(i)) + peak_interval_width/2, max(freq));
        peak_intervals(:, i) = [lower_bound, upper_bound];
    end
end

