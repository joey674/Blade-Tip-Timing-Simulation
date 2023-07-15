function Magn_smooth = damping_re_smooth(Magn,locs)
    window_size = 100;   
    Magn_smooth = Magn;
    for i = 1:length(Magn)
        if ~ismember(i,locs)
            window_start = max(1,i-window_size);
            window_end = min(length(Magn),i+window_size);
            Magn_smooth(i) = mean(Magn(window_start:window_end));
        end 
    end    
end