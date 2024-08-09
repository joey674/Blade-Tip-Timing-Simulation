function [fig, peaks_idx] = FindPeakManual(fig, freq, magn, peaks_idx)
    figure(fig); 
    hold on;

    %{
        plot exist peaks
    %} 
    for i = 1:length(peaks_idx)
        plot(freq(peaks_idx(i)), magn(peaks_idx(i)), 'bo', 'DisplayName', ['Peak idx ' num2str(peaks_idx(i))]);
    end
    
    %{
        add peaks manually
    %} 
    while true
        [freq_get, ~] = ginput(1); 
        if isempty(freq_get)  
            break;
        end
        [~, idx] = min(abs(freq - freq_get));
        idx = round(idx);
        searchRange = 30;
        startIndex = max(1, idx - searchRange);
        endIndex = min(length(magn), idx + searchRange);
        [~, maxIndex] = max(magn(startIndex:endIndex));
        idx = startIndex + maxIndex - 1;
        plot(freq(idx), magn(idx), 'bo', 'DisplayName', ['Peak idx ' num2str(idx)]);
        peaks_idx = [peaks_idx, idx];            
    end
end


