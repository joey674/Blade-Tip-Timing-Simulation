function n_modes = MDOF_FindNumOfModes(blade,smoothness)
peak_prominence = 5;
peak_interval_width = 10;

n_blades = length(blade);
for i = 1:n_blades
    lens(i) = length([blade{i}.magn]);
end
magn_all = NaN(max(lens),n_blades);
for i = 1:n_blades
    magn_all(1:length([blade{i}.magn]),i) = [blade{i}.magn];
end
magn_max = max(magn_all,[],2,"omitnan");
magn_max = smoothdata(magn_max,"movmean",smoothness);

peak_prominence = 0.2 * (max(magn_max) - mean(magn_max));

% find peak interval
[pks, locs] = findpeaks(magn_max, 'MinPeakProminence', peak_prominence); 
n_modes = length(pks);
peak_intervals = zeros(n_modes,2);
for i = 1:n_modes
    lower_bound = blade{1}(locs(i)).freq - peak_interval_width/2;
    upper_bound = blade{1}(locs(i)).freq + peak_interval_width/2;
    peak_intervals(i, :) = [lower_bound, upper_bound];
end

figure('units','normalized','outerposition',[0 0 0.7 0.7]);  
plot([blade{1}.freq],magn_max,'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Magnitude');
hold on;
plot([blade{1}(locs).freq], pks, 'ro');
hold off;
fprintf("prominence:%d,interval_width:%d,n_modes:%d\n",peak_prominence,peak_interval_width,n_modes);
end



