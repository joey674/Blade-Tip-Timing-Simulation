function magn_smoothed = ReduceNoise(magn, freq, err)

    magn_downsampled = smoothdata(magn, 'sgolay', 200);

    magn_smoothed = magn_downsampled;  
end


