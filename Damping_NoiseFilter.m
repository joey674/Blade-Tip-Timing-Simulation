function filtered_data = Damping_NoiseFilter(data)   
    data_movmean = smoothdata(data,'movmean',100);
    % data_gauss = smoothdata(data,'gaussian',500);
    % conf = max(data_movmean)/max(data_gauss);
    % filtered_data = data_gauss*conf;
    filtered_data = data_movmean;
end
