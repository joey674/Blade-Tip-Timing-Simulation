function uniform_blade = Interpolate(blade)
     for blade_idx = 1:12
        %% Initialize
        blade_data = blade{blade_idx};
        freq = [blade_data.freq];
        magn = [blade_data.magn];
        phase = [blade_data.phase];
        err = [blade_data.err]; 

        %% Interpolate to the same distance
        min_freq = min(freq);
        max_freq = max(freq);
        N = length(freq); 
        
        uniform_freq = linspace(min_freq, max_freq, N);
        
        % interpolate magn, phase, err 
        uniform_magn = interp1(freq, magn, uniform_freq, 'linear'); 
        uniform_phase = interp1(freq, phase, uniform_freq, 'linear');
        uniform_err = interp1(freq, err, uniform_freq, 'linear');

        uniform_blade_data.freq = uniform_freq;
        uniform_blade_data.magn = uniform_magn;
        uniform_blade_data.phase = uniform_phase;
        uniform_blade_data.err = uniform_err;

        uniform_blade{blade_idx} = uniform_blade_data;
     end
end

